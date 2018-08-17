#----------------------------------------------------------------------------------
# In order to pass the CRAN submission process, we currently cannot put h2oEnsemble 
# in the required/ suggested list of libraries. But once it's available in R, we 
# will add h2oEnsemble back into the list.
# For details, please visit https://www.stat.berkeley.edu/~ledell/R/h2oEnsemble.pdf
# https://github.com/h2oai/h2o-3/tree/master/h2o-r/ensemble/h2oEnsemble-package
# In this R file, all required functions from h2oEnsemble package has been included.
#----------------------------------------------------------------------------------

#--------------------------------------- h2o.ensemble --------------------------------------
# Purpose: This function creates a "Super Learner" (stacking) ensemble using the 
# H2O base learning algorithms specified by the user
#-------------------------------------------------------------------------------------------
h2o.ensemble <- function(x, y, training_frame, 
                         model_id = NULL, validation_frame = NULL,
                         family = c("AUTO", "binomial", "gaussian", "quasibinomial", "poisson", "gamma", "tweedie", "laplace", "quantile", "huber"),
                         learner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper", "h2o.gbm.wrapper", "h2o.deeplearning.wrapper"),
                         metalearner = "h2o.glm.wrapper",
                         cvControl = list(V = 5, shuffle = TRUE),  #maybe change this to cv_control
                         seed = 1,
                         parallel = "seq",  #only seq implemented
                         keep_levelone_data = TRUE) 
{
  
  starttime <- Sys.time()
  runtime <- list()
  
  # Training_frame may be a key or an H2O H2OFrame object
  if ((!inherits(training_frame, "Frame") && !inherits(training_frame, "H2OFrame")))
    tryCatch(training_frame <- h2o.getFrame(training_frame),
             error = function(err) {
               stop("argument \"training_frame\" must be a valid H2OFrame or id")
             })
  if (!is.null(validation_frame)) {
    if (is.character(validation_frame))
      tryCatch(validation_frame <- h2o.getFrame(validation_frame),
               error = function(err) {
                 stop("argument \"validation_frame\" must be a valid H2OFrame or id")
               })
  }
  N <- dim(training_frame)[1L]  #Number of observations in training set
  if (is.null(validation_frame)) {
    validation_frame <- training_frame
  }
  
  # Determine prediction task family type automatically
  # TO DO: Add auto-detection for other distributions like gamma - right now auto-detect as "gaussian"
  family <- match.arg(family)
  if (family == "AUTO") {
    if (is.factor(training_frame[,y])) {
      numcats <- length(h2o.levels(training_frame[,y]))
      if (numcats == 2) {
        family <- "binomial" 
      } else {
        stop("Multinomial case not yet implemented for h2o.ensemble. Check here for progress: https://0xdata.atlassian.net/browse/PUBDEV-2355")
      }
    } else {
      family <- "gaussian"
    }
  }
  # Check that if specified, family matches data type for response
  # binomial must be factor/enum and gaussian must be numeric
  if (family %in% c("gaussian", "quasibinomial", "poisson", "gamma", "tweedie", "laplace", "quantile", "huber")) {
    if (!is.numeric(training_frame[,y])) {
      stop("When `family` is one of {gaussian, quasibinomial, poisson, gamma, tweedie, laplace, quantile, huber}, the repsonse column must be numeric.")
    }
    # TO DO: Update this ylim calc when h2o.range method gets implemented for H2OFrame cols
    #ylim <- c(min(training_frame[,y]), max(training_frame[,y]))  #Used to enforce bounds
    ylim <- h2o.range(training_frame[,y])
    if (family == "gamma") {
      if (ylim[1] <= 0) {
        stop("family = gamma requires a positive respone")
      }      
    }
  } else {
    if (!is.factor(training_frame[,y])) {
      stop("When `family` is binomial, the repsonse column must be a factor.")
    } else {
      numcats <- length(h2o.levels(training_frame[,y]))
      if (numcats > 2) {
        stop("Multinomial case not yet implemented for h2o.ensemble. Check here for progress: https://0xdata.atlassian.net/browse/PUBDEV-2355")
      } 
    }
    ylim <- NULL
  }
  
  # Update control args by filling in missing list elements
  cvControl <- do.call(".cv_control", cvControl)
  V <- cvControl$V      #Number of CV folds
  L <- length(learner)  #Number of distinct learners
  idxs <- expand.grid(1:V,1:L)
  names(idxs) <- c("v","l")


  # Validate learner and metalearner arguments
  if (length(metalearner)>1 | !is.character(metalearner) | !exists(metalearner)) {
    stop("The 'metalearner' argument must be a string, specifying the name of a base learner wrapper function.")
  }
  if (sum(!sapply(learner, exists))>0) {
    stop("'learner' function name(s) not found.")
  }
  if (!exists(metalearner)) {
    stop("'metalearner' function name not found.")
  }
    
  # Check remaining args
  if (inherits(parallel, "character")) {
    if (!(parallel %in% c("seq","multicore"))) {
      stop("'parallel' must be either 'seq' or 'multicore' or a snow cluster object")
    }
  } else if (!inherits(parallel, "cluster")) {
    stop("'parallel' must be either 'seq' or 'multicore' or a snow cluster object")
  }
  
  # Begin ensemble code
  if (is.numeric(seed)) set.seed(seed)  #If seed is specified, set seed prior to next step
  folds <- sample(rep(seq(V), ceiling(N/V)))[1:N]  # Cross-validation folds (stratified folds not yet supported)
  training_frame$fold_id <- as.h2o(folds)  # Add a fold_id column for each observation so we can subset by row later

  # What type of metalearning function do we have?
  # The h2o version is memory-optimized (the N x L level-one matrix, Z, never leaves H2O memory);
  # SuperLearner metalearners provide additional metalearning algos, but has a much bigger memory footprint
  if (grepl("^SL.", metalearner)) {
    metalearner_type <- "SuperLearner"
  } else if (grepl("^h2o.", metalearner)){
    metalearner_type <- "h2o"
  }
  
  # Create the Z matrix of cross-validated predictions
  mm <- .make_Z(x = x, y = y, training_frame = training_frame, 
                family = family, 
                learner = learner, 
                parallel = parallel, 
                seed = seed, 
                V = V, 
                L = L, 
                idxs = idxs,
                metalearner_type = metalearner_type)
  # TO DO: Could pass on the metalearner arg instead of metalearner_type and get this info internally
  basefits <- mm$basefits
  Z <- mm$Z  #pure Z (dimension N x L)
  
  # Metalearning: Regress y onto Z to learn optimal combination of base models
  # TO DO: Replace grepl for metalearner_type
  # TO DO: Pass on additional args to match.fun(metalearner) for h2o type
  print("Metalearning")
  if (is.numeric(seed)) set.seed(seed)  #If seed given, set seed prior to next step
  if (grepl("^SL.", metalearner)) {
    # this is very hacky and should be used only for testing
    if (is.character(family)) {
      familyFun <- get(family, mode = "function", envir = parent.frame())
      #print(familyFun$family)  #does not work for SL.glmnet
    } 
    Zdf <- as.data.frame(Z)
    Y <- as.data.frame(training_frame[,c(y)])[,1]
    # TO DO: for parity, need to add y col to Z like we do below
    runtime$metalearning <- system.time(metafit <- match.fun(metalearner)(Y = Y, 
                                                                          X = Zdf, 
                                                                          newX = Zdf, 
                                                                          family = familyFun, 
                                                                          id = seq(N), 
                                                                          obsWeights = rep(1,N)), gcFirst = FALSE)
  } else {
    Z$y <- training_frame[,c(y)]  # do we want to add y to the Z frame?  
    runtime$metalearning <- system.time(metafit <- match.fun(metalearner)(x = learner, 
                                                                          y = "y", 
                                                                          training_frame = Z, 
                                                                          validation_frame = NULL, 
                                                                          family = family), gcFirst = FALSE)
  }
  
  # Since baselearning is now performed along with CV, see if we can get this info, or deprecate this
  runtime$baselearning <- NULL
  runtime$total <- Sys.time() - starttime
  
  # Keep level-one data?
  if (!keep_levelone_data) {
    Z <- NULL
  }
  
  # Ensemble model
  out <- list(x = x,
              y = y, 
              family = family, 
              learner = learner,
              metalearner = metalearner,
              cvControl = cvControl,
              folds = folds,
              ylim = ylim, 
              seed = seed,
              parallel = parallel,
              basefits = basefits, 
              metafit = metafit,
              levelone = Z,  #levelone = cbind(Z, y)
              runtime = runtime,
              h2o_version = packageVersion(pkg = "h2o"),
              h2oEnsemble_version = packageVersion(pkg = "h2oEnsemble"))
  class(out) <- "h2o.ensemble"
  return(out)
}


#----------------------------------------- .make_Z -----------------------------------------
# Purpose: Generate the CV predicted values for all learners
#-------------------------------------------------------------------------------------------
.make_Z <- function(x, y, training_frame, family, learner, parallel, seed, V, L, idxs, metalearner_type = c("h2o", "SuperLearner")) {
  
  # Wrapper function for .fitFun to record system.time
  .fitWrapper <- function(l, y, xcols, training_frame, validation_frame, family, learner, seed, fold_column) {
    print(sprintf("Cross-validating and training base learner %s: %s", l, learner[l]))
    fittime <- system.time(fit <- .fitFun(l, y, xcols, training_frame, validation_frame, family, 
                                          learner, seed, fold_column), gcFirst=FALSE)
    return(list(fit=fit, fittime=fittime))
  }
  
  # Do V-fold cross-validation of each learner (in a loop/apply over 1:L)...
  fitlist <- sapply(X = 1:L, FUN = .fitWrapper, y = y, xcols = x, training_frame = training_frame,
                    validation_frame = NULL, family = family, learner = learner, 
                    seed = seed, fold_column = "fold_id", 
                    simplify = FALSE)
  
  runtime <- list()
  runtime$cv <- lapply(fitlist, function(ll) ll$fittime)
  names(runtime$cv) <- learner
  basefits <- lapply(fitlist, function(ll) ll$fit)  #Base fits (trained on full data) to be saved
  names(basefits) <- learner      
  
  # In the case of binary classification, a 3-col HDF is returned, colnames == c("predict", "p0", "p1")
  # In the case of regression, 1-col HDF is already returned, colname == "predict"
  .compress_cvpred_into_1col <- function(l, family) {
    # return the frame_id of the resulting 1-col Hdf of cvpreds for learner l
    if (family %in% c("bernoulli", "binomial")) {
      predlist <- sapply(1:V, function(v) h2o.getFrame(basefits[[l]]@model$cross_validation_predictions[[v]]$name)[,3], simplify = FALSE)
    } else {
      predlist <- sapply(1:V, function(v) h2o.getFrame(basefits[[l]]@model$cross_validation_predictions[[v]]$name)$predict, simplify = FALSE)
    }
    cvpred_sparse <- h2o.cbind(predlist)  #N x V Hdf with rows that are all zeros, except corresponding to the v^th fold if that rows is associated with v
    cvpred_col <- apply(cvpred_sparse, 1, sum)
    return(cvpred_col)
  } 
  cvpred_framelist <- sapply(1:L, function(l) .compress_cvpred_into_1col(l, family))
  Z <- h2o.cbind(cvpred_framelist)
  names(Z) <- learner
  return(list(Z = Z, basefits = basefits))
}


#----------------------------------------- .fitFun -----------------------------------------
# Purpose: Train a model using learner l 
#-------------------------------------------------------------------------------------------
.fitFun <- function(l, y, x, training_frame, validation_frame, family, learner, seed, fold_column) {
  if (!is.null(fold_column)) cv = TRUE
  if (is.numeric(seed)) set.seed(seed)  #If seed given, set seed prior to next step
  if (("x" %in% names(formals(learner[l]))) && (as.character(formals(learner[l])$x)[1] != "")) {
    # Special case where we pass a subset of the colnames, x, in a custom learner function wrapper
    fit <- match.fun(learner[l])(y = y, training_frame = training_frame, validation_frame = NULL, family = family, fold_column = fold_column, keep_cross_validation_folds = TRUE)
  } else {
    # Use all predictors in training set for training
    fit <- match.fun(learner[l])(y = y, x = x, training_frame = training_frame, validation_frame = NULL, family = family, fold_column = fold_column, keep_cross_validation_folds = TRUE)
  }
  #fit <- get(learner[l], mode = "function", envir = parent.frame())(y = y, x = x, training_frame = training_frame, validation_frame = NULL, family = family, fold_column = fold_column, keep_cross_validation_folds = cv)
  return(fit)
}

                             
#---------------------------------------- .cv_control --------------------------------------
# Purpose: Parameters that control the CV process
#-------------------------------------------------------------------------------------------  
.cv_control <- function(V = 5L, stratifyCV = TRUE, shuffle = TRUE){
  # Parameters that control the CV process
  # Only part of this being used currently --  
  # Stratification is not enabled yet in the h2o.ensemble function.
  # We can use a modified SuperLearner::CVFolds function (or similar) to 
  # enable stratification by outcome in the future.
  
  V <- as.integer(V)  #Number of cross-validation folds
  if(!is.logical(stratifyCV)) {
    stop("'stratifyCV' must be logical")
  }
  if(!is.logical(shuffle)) {
    stop("'shuffle' must be logical")
  }  
  return(list(V = V, stratifyCV = stratifyCV, shuffle = shuffle))
}


#----------------------------------- predict.h2o.ensemble ----------------------------------
# Purpose: Predict method for an ’h2o.ensemble’ object
#-------------------------------------------------------------------------------------------                             
predict.h2o.ensemble <- function(object, newdata, ...) {
  
  if (object$family == "binomial") {
    basepred <- h2o.cbind(sapply(object$basefits, function(ll) h2o.predict(object = ll, newdata = newdata)[,3]))
  } else {
    basepred <- h2o.cbind(sapply(object$basefits, function(ll) h2o.predict(object = ll, newdata = newdata)[,1]))
  }
  names(basepred) <- names(object$basefits)
  
  if (grepl("H2O", class(object$metafit))) {
    # H2O ensemble metalearner from wrappers.R
    pred <- h2o.predict(object = object$metafit, newdata = basepred)
  } else {
    # SuperLearner wrapper function metalearner
    basepreddf <- as.data.frame(basepred)  
    pred <- predict(object = object$metafit$fit, newdata = basepreddf)
  }
  out <- list(pred = pred, basepred = basepred)
  return(out)
}

                                 
#-------------------------------------------------------------------------------------------                                                              
# Set of default wrappers to create a uniform interface for h2o supervised ML functions (H2O 3.0 and above)
# These wrapper functions should always be compatible with the master branch of: https://github.com/h2oai/h2o-3

# Example of a wrapper function:
h2o.example.wrapper <- function(x, y, training_frame, model_id = NULL, family = c("gaussian", "binomial"), ...) {
  # This function is just an example.  
  # You can wrap any H2O learner inside a wrapper function, example: h2o.glm
  h2o.glm(x = x, y = y, training_frame = training_frame, family = family)
}

# H2O Algorithm function wrappers for:
# h2o.glm
# h2o.randomForest
# h2o.gbm
# h2o.deeplearning

h2o.glm.wrapper <- function(x, y, training_frame, 
                            model_id = NULL,
                            validation_frame = NULL,
                            nfolds = 0,
                            seed = -1,
                            keep_cross_validation_predictions = TRUE,
                            keep_cross_validation_fold_assignment = FALSE,
                            fold_assignment = c("AUTO", "Random", "Modulo", "Stratified"),
                            fold_column = NULL,
                            ignore_const_cols = TRUE,
                            score_each_iteration = FALSE,
                            offset_column = NULL,
                            weights_column = NULL,
                            family = "AUTO",
                            #family = c("AUTO", "binomial", "gaussian", "quasibinomial", "poisson", "gamma", "tweedie", "laplace", "quantile", "huber"),
                            tweedie_variance_power = 0,
                            tweedie_link_power = 1,
                            solver = c("AUTO", "IRLSM", "L_BFGS", "COORDINATE_DESCENT_NAIVE", "COORDINATE_DESCENT"),
                            alpha = NULL,
                            lambda = NULL,
                            lambda_search = FALSE,
                            early_stopping = TRUE,
                            nlambdas = -1,
                            standardize = TRUE,
                            missing_values_handling = c("MeanImputation", "Skip"),
                            compute_p_values = FALSE,
                            remove_collinear_columns = FALSE,
                            intercept = TRUE,
                            non_negative = FALSE,
                            max_iterations = -1,
                            objective_epsilon = -1,
                            beta_epsilon = 0.0001,
                            gradient_epsilon = -1,
                            link = c("family_default", "identity", "logit", "log", "inverse", "tweedie"),
                            prior = -1,
                            lambda_min_ratio = -1,
                            beta_constraints = NULL,
                            max_active_predictors = -1,
                            interactions = NULL,
                            balance_classes = FALSE,
                            class_sampling_factors = NULL,
                            max_after_balance_size = 5.0,
                            max_hit_ratio_k = 0,
                            max_runtime_secs = 0, ...) {
  
  # If family is not specified, set it using the datatype of the response column
  #family <- match.arg(family)
  if (family == "AUTO") {
    if (is.factor(training_frame[,y])) {
      family <- "binomial"
    } else {
      family <- "gaussian"
    }
  } else if (family %in% c("laplace", "quantile", "huber")) { # not supported by GLM
      family <- "gaussian"
  }
  
  # Also, offset_column, weights_column, intercept not implemented at the moment due to similar bug as beta_constraints
  # intercept argument not currently supported due to GLM bug with explicitly setting interactions = NULL (the default) 
  h2o.glm(x = x, 
          y = y, 
          training_frame = training_frame, 
          model_id = model_id, 
          validation_frame = validation_frame,
          nfolds = nfolds,
          seed = seed,
          keep_cross_validation_predictions = TRUE,  #must have for stacking
          keep_cross_validation_fold_assignment = keep_cross_validation_fold_assignment, 
          fold_assignment = match.arg(fold_assignment),
          fold_column = fold_column,
          ignore_const_cols = ignore_const_cols,
          score_each_iteration = score_each_iteration,
          offset_column = offset_column,
          weights_column = weights_column,
          family = family, 
          tweedie_variance_power = tweedie_variance_power,
          tweedie_link_power = tweedie_link_power,        
          solver = match.arg(solver), 
          alpha = alpha, 
          lambda = lambda, 
          lambda_search = lambda_search,
          early_stopping = early_stopping,
          nlambdas = nlambdas,           
          standardize = standardize,
          missing_values_handling = match.arg(missing_values_handling),
          compute_p_values = compute_p_values,
          remove_collinear_columns = remove_collinear_columns,
          intercept = intercept,
          non_negative = non_negative, 
          max_iterations = max_iterations,
          objective_epsilon = objective_epsilon,
          beta_epsilon = beta_epsilon, 
          gradient_epsilon = gradient_epsilon, 
          link = match.arg(link),
          prior = prior,
          lambda_min_ratio = lambda_min_ratio,
          beta_constraints = beta_constraints,
          max_active_predictors = max_active_predictors,
          #interactions = interactions,  #causes a bug when set to NULL (the default), the h2o.glm function needs to be fixed: https://0xdata.atlassian.net/browse/PUBDEV-4698
          balance_classes = balance_classes,
          class_sampling_factors = class_sampling_factors,
          max_after_balance_size = max_after_balance_size,
          max_hit_ratio_k = max_hit_ratio_k, 
          max_runtime_secs = max_runtime_secs)
}



h2o.gbm.wrapper <- function(x, y, training_frame, model_id = NULL,
                            family = "AUTO",
                            validation_frame = NULL,
                            nfolds = 0,
                            keep_cross_validation_predictions = FALSE,
                            keep_cross_validation_fold_assignment = FALSE,
                            score_each_iteration = FALSE,
                            score_tree_interval = 0,
                            fold_assignment = c("AUTO", "Random", "Modulo", "Stratified"),
                            fold_column = NULL,
                            ignore_const_cols = TRUE,
                            offset_column = NULL,
                            weights_column = NULL,
                            balance_classes = FALSE,
                            class_sampling_factors = NULL,
                            max_after_balance_size = 5.0,
                            max_hit_ratio_k = 0,
                            ntrees = 50,
                            max_depth = 5,
                            min_rows = 10,
                            nbins = 20,
                            nbins_top_level = 1024,
                            nbins_cats = 1024,
                            #r2_stopping = 1.797693135e+308,  #deprecated
                            stopping_rounds = 0,
                            stopping_metric = c("AUTO", "deviance", "logloss", "MSE", "RMSE", "MAE", "RMSLE", "AUC", "lift_top_group", "misclassification", "mean_per_class_error"),
                            stopping_tolerance = 0.001,
                            max_runtime_secs = 0,
                            seed = -1,
                            build_tree_one_node = FALSE,
                            learn_rate = 0.1,
                            learn_rate_annealing = 1,
                            #distribution = c("AUTO", "bernoulli", "multinomial", "gaussian", "poisson", "gamma", "tweedie", "laplace", "quantile", "huber"),
                            quantile_alpha = 0.5,
                            tweedie_power = 1.5,
                            huber_alpha = 0.9,
                            checkpoint = NULL,
                            sample_rate = 1,
                            sample_rate_per_class = NULL,
                            col_sample_rate = 1,
                            col_sample_rate_change_per_level = 1,
                            col_sample_rate_per_tree = 1,
                            min_split_improvement = 1e-05,
                            histogram_type = c("AUTO", "UniformAdaptive", "Random", "QuantilesGlobal", "RoundRobin"),
                            max_abs_leafnode_pred = 1.797693135e+308,
                            pred_noise_bandwidth = 0,
                            categorical_encoding = c("AUTO", "Enum", "OneHotInternal", "OneHotExplicit", "Binary", "Eigen", "LabelEncoder", "SortByResponse", "EnumLimited"),
                            calibrate_model = FALSE,
                            calibration_frame = NULL,
                            #verbose = FALSE,  #remove so that this is compatible with earlier versions of H2O
                            ...) {
  
  # If family is not specified, set it using the datatype of the response column
  # GBM uses `distribution` instead of `family` so we set `distribution` here instead
  if (family == "AUTO") {
    if (is.factor(training_frame[,y])) {
      distribution <- "bernoulli"
    } else {
      distribution <- "gaussian"
    }
  } else if (family %in% c("laplace", "quantile", "huber")) { # extra distributions for GBM
    distribution <- family
  } else if (family == "binomial") {
    distribution <- "bernoulli"
  } else {
    # not supported by GBM, so we set to "gaussian"
    distribution <- "gaussian"
  }
  
  h2o.gbm(x = x, 
          y = y, 
          training_frame = training_frame, 
          model_id = model_id,
          validation_frame = validation_frame,
          nfolds = nfolds,
          keep_cross_validation_predictions = TRUE,  #must have for stacking
          keep_cross_validation_fold_assignment = keep_cross_validation_fold_assignment,
          score_each_iteration = score_each_iteration,
          score_tree_interval = score_tree_interval,
          fold_assignment = match.arg(fold_assignment),
          fold_column = fold_column,
          ignore_const_cols = ignore_const_cols,
          offset_column = offset_column,
          weights_column = weights_column,
          balance_classes = balance_classes,
          class_sampling_factors = class_sampling_factors,
          max_after_balance_size = max_after_balance_size,
          max_hit_ratio_k = max_hit_ratio_k,
          ntrees = ntrees,
          max_depth = max_depth,
          min_rows = min_rows,
          nbins = nbins,
          nbins_top_level = nbins_top_level,
          nbins_cats = nbins_cats,
          #r2_stopping = r2_stopping,
          stopping_rounds = stopping_rounds,
          stopping_metric = match.arg(stopping_metric),
          stopping_tolerance = stopping_tolerance,
          max_runtime_secs = max_runtime_secs,
          seed = seed,
          build_tree_one_node = build_tree_one_node,
          learn_rate = learn_rate,
          learn_rate_annealing = learn_rate_annealing,
          distribution = distribution, # set above
          quantile_alpha = quantile_alpha,
          tweedie_power = tweedie_power,
          huber_alpha = huber_alpha,
          checkpoint = checkpoint,
          sample_rate = sample_rate,
          sample_rate_per_class = sample_rate_per_class,
          col_sample_rate = col_sample_rate,
          col_sample_rate_change_per_level = col_sample_rate_change_per_level,
          col_sample_rate_per_tree = col_sample_rate_per_tree,
          min_split_improvement = min_split_improvement,
          histogram_type = match.arg(histogram_type),
          max_abs_leafnode_pred = max_abs_leafnode_pred,
          pred_noise_bandwidth = pred_noise_bandwidth,
          categorical_encoding = match.arg(categorical_encoding),
          calibrate_model = calibrate_model,
          calibration_frame = calibration_frame)
}


h2o.randomForest.wrapper <- function(x, y, training_frame, model_id = NULL,
                                     validation_frame = NULL,
                                     nfolds = 0,
                                     keep_cross_validation_predictions = TRUE,
                                     keep_cross_validation_fold_assignment = FALSE,
                                     score_each_iteration = FALSE,
                                     score_tree_interval = 0,
                                     fold_assignment = c("AUTO", "Random", "Modulo", "Stratified"),
                                     fold_column = NULL,
                                     ignore_const_cols = TRUE,
                                     offset_column = NULL,
                                     weights_column = NULL,
                                     balance_classes = FALSE,
                                     class_sampling_factors = NULL,
                                     max_after_balance_size = 5.0,
                                     max_hit_ratio_k = 0,
                                     ntrees = 50,
                                     max_depth = 20,
                                     min_rows = 1,
                                     nbins = 20,
                                     nbins_top_level = 1024,
                                     nbins_cats = 1024,
                                     #r2_stopping = 1.797693135e+308,  #deprecated
                                     stopping_rounds = 0,
                                     stopping_metric = c("AUTO", "deviance", "logloss", "MSE", "RMSE", "MAE", "RMSLE", "AUC", "lift_top_group", "misclassification", "mean_per_class_error"),
                                     stopping_tolerance = 0.001,
                                     max_runtime_secs = 0,
                                     seed = -1,
                                     build_tree_one_node = FALSE,
                                     mtries = -1,
                                     sample_rate = 0.6320000291,
                                     sample_rate_per_class = NULL,
                                     binomial_double_trees = FALSE,
                                     checkpoint = NULL,
                                     col_sample_rate_change_per_level = 1,
                                     col_sample_rate_per_tree = 1,
                                     min_split_improvement = 1e-05,
                                     histogram_type = c("AUTO", "UniformAdaptive", "Random", "QuantilesGlobal", "RoundRobin"),
                                     categorical_encoding = c("AUTO", "Enum", "OneHotInternal", "OneHotExplicit", "Binary", "Eigen", "LabelEncoder", "SortByResponse", "EnumLimited"),
                                     calibrate_model = FALSE,
                                     calibration_frame = NULL,
                                     #verbose = FALSE,  #remove so that this is compatible with earlier versions of H2O
                                     ...) {
  
  # Currently ignoring the `family` arg (so it's removed), will get class from outcome in H2OFrame
  # TO DO: Add a check to make sure that outcome/family type is consistent with specified family
  h2o.randomForest(x = x, 
                   y = y, 
                   training_frame = training_frame, 
                   model_id = model_id, 
                   validation_frame = validation_frame,
                   nfolds = nfolds,
                   keep_cross_validation_predictions = TRUE,  #must have for stacking
                   keep_cross_validation_fold_assignment = keep_cross_validation_fold_assignment,
                   score_each_iteration = score_each_iteration,
                   score_tree_interval = score_tree_interval,
                   fold_assignment = match.arg(fold_assignment),
                   fold_column = fold_column,
                   ignore_const_cols = ignore_const_cols,
                   offset_column = offset_column,
                   weights_column = weights_column,
                   balance_classes = balance_classes,
                   class_sampling_factors = class_sampling_factors,
                   max_after_balance_size = max_after_balance_size,
                   max_hit_ratio_k = max_hit_ratio_k,
                   ntrees = ntrees,
                   max_depth = max_depth,
                   min_rows = min_rows,
                   nbins = nbins,
                   nbins_top_level = nbins_top_level,
                   nbins_cats = nbins_cats,
                   #r2_stopping = r2_stopping,
                   stopping_rounds = stopping_rounds,
                   stopping_metric = match.arg(stopping_metric),
                   stopping_tolerance = stopping_tolerance,
                   max_runtime_secs = max_runtime_secs,
                   seed = seed,
                   build_tree_one_node = build_tree_one_node,
                   mtries = mtries,
                   sample_rate = sample_rate,
                   sample_rate_per_class = sample_rate_per_class,
                   binomial_double_trees = binomial_double_trees,
                   checkpoint = checkpoint,
                   col_sample_rate_change_per_level = col_sample_rate_change_per_level,
                   col_sample_rate_per_tree = col_sample_rate_per_tree,
                   min_split_improvement = min_split_improvement,
                   histogram_type = match.arg(histogram_type),
                   categorical_encoding = match.arg(categorical_encoding),
                   calibrate_model = calibrate_model,
                   calibration_frame = calibration_frame)
}


h2o.deeplearning.wrapper <- function(x, y, training_frame, model_id = NULL,
                                     family = "AUTO", 
                                     validation_frame = NULL,
                                     nfolds = 0,
                                     keep_cross_validation_predictions = FALSE,
                                     keep_cross_validation_fold_assignment = FALSE,
                                     fold_assignment = c("AUTO", "Random", "Modulo", "Stratified"),
                                     fold_column = NULL,
                                     ignore_const_cols = TRUE,
                                     score_each_iteration = FALSE,
                                     weights_column = NULL,
                                     offset_column = NULL,
                                     balance_classes = FALSE,
                                     class_sampling_factors = NULL,
                                     max_after_balance_size = 5.0,
                                     max_hit_ratio_k = 0,
                                     checkpoint = NULL,
                                     pretrained_autoencoder = NULL,
                                     overwrite_with_best_model = TRUE,
                                     use_all_factor_levels = TRUE,
                                     standardize = TRUE,
                                     activation = c("Tanh", "TanhWithDropout", "Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout"),
                                     hidden = c(200, 200),
                                     epochs = 10,
                                     train_samples_per_iteration = -2,
                                     target_ratio_comm_to_comp = 0.05,
                                     seed = -1,
                                     adaptive_rate = TRUE,
                                     rho = 0.99,
                                     epsilon = 1e-08,
                                     rate = 0.005,
                                     rate_annealing = 1e-06,
                                     rate_decay = 1,
                                     momentum_start = 0,
                                     momentum_ramp = 1000000,
                                     momentum_stable = 0,
                                     nesterov_accelerated_gradient = TRUE,
                                     input_dropout_ratio = 0,
                                     hidden_dropout_ratios = NULL,
                                     l1 = 0,
                                     l2 = 0,
                                     max_w2 = 3.4028235e+38,
                                     initial_weight_distribution = c("UniformAdaptive", "Uniform", "Normal"),
                                     initial_weight_scale = 1,
                                     initial_weights = NULL,
                                     initial_biases = NULL,
                                     loss = c("Automatic", "CrossEntropy", "Quadratic", "Huber", "Absolute", "Quantile"),
                                     #distribution = c("AUTO", "bernoulli", "multinomial", "gaussian", "poisson", "gamma", "tweedie", "laplace", "quantile", "huber"),
                                     quantile_alpha = 0.5,
                                     tweedie_power = 1.5,
                                     huber_alpha = 0.9,
                                     score_interval = 5,
                                     score_training_samples = 10000,
                                     score_validation_samples = 0,
                                     score_duty_cycle = 0.1,
                                     classification_stop = 0,
                                     regression_stop = 1e-06,
                                     stopping_rounds = 5,
                                     stopping_metric = c("AUTO", "deviance", "logloss", "MSE", "RMSE", "MAE", "RMSLE", "AUC", "lift_top_group", "misclassification", "mean_per_class_error"),
                                     stopping_tolerance = 0,
                                     max_runtime_secs = 0,
                                     score_validation_sampling = c("Uniform", "Stratified"),
                                     diagnostics = TRUE,
                                     fast_mode = TRUE,
                                     force_load_balance = TRUE,
                                     variable_importances = TRUE,
                                     replicate_training_data = TRUE,
                                     single_node_mode = FALSE,
                                     shuffle_training_data = FALSE,
                                     missing_values_handling = c("MeanImputation", "Skip"),
                                     quiet_mode = FALSE,
                                     autoencoder = FALSE,
                                     sparse = FALSE,
                                     col_major = FALSE,
                                     average_activation = 0,
                                     sparsity_beta = 0,
                                     max_categorical_features = 2147483647,
                                     reproducible = FALSE,
                                     export_weights_and_biases = FALSE,
                                     mini_batch_size = 1,
                                     categorical_encoding = c("AUTO", "Enum", "OneHotInternal", "OneHotExplicit", "Binary", "Eigen", "LabelEncoder", "SortByResponse", "EnumLimited"),
                                     elastic_averaging = FALSE,
                                     elastic_averaging_moving_rate = 0.9,
                                     elastic_averaging_regularization = 0.001,
                                     #verbose = FALSE,  #remove so that this is compatible with earlier versions of H2O
                                     ...) {
  
  # If family is not specified, set it using the datatype of the response column
  # GBM uses `distribution` instead of `family` so we set `distribution` here instead
  if (family == "AUTO") {
    if (is.factor(training_frame[,y])) {
      distribution <- "bernoulli"
    } else {
      distribution <- "gaussian"
    }
  } else if (family %in% c("laplace", "quantile", "huber")) { # extra distributions for DL
    distribution <- family
  } else if (family == "binomial") {
    distribution <- "bernoulli"
  } else {
    # not supported by DL, so we set to "gaussian"
    distribution <- "gaussian"
  }
  
  h2o.deeplearning(x = x, 
                   y = y, 
                   training_frame = training_frame, 
                   model_id = model_id,
                   validation_frame = NULL,
                   nfolds = nfolds,
                   keep_cross_validation_predictions = TRUE,  #must have for stacking
                   keep_cross_validation_fold_assignment = keep_cross_validation_fold_assignment,
                   fold_assignment = match.arg(fold_assignment),
                   fold_column = fold_column,
                   ignore_const_cols = ignore_const_cols,
                   score_each_iteration = score_each_iteration,
                   weights_column = weights_column,
                   offset_column = offset_column,
                   balance_classes = balance_classes,
                   class_sampling_factors = class_sampling_factors,
                   max_after_balance_size = max_after_balance_size,
                   max_hit_ratio_k = max_hit_ratio_k,
                   checkpoint = checkpoint,
                   pretrained_autoencoder = pretrained_autoencoder,
                   overwrite_with_best_model = overwrite_with_best_model,
                   use_all_factor_levels = use_all_factor_levels,
                   standardize = standardize,
                   activation = match.arg(activation),
                   hidden = hidden,
                   epochs = epochs,
                   train_samples_per_iteration = train_samples_per_iteration,
                   target_ratio_comm_to_comp = target_ratio_comm_to_comp,
                   seed = seed,
                   adaptive_rate = adaptive_rate,
                   rho = rho,
                   epsilon = epsilon,
                   rate = rate,
                   rate_annealing = rate_annealing,
                   rate_decay = rate_decay,
                   momentum_start = momentum_start,
                   momentum_ramp = momentum_ramp,
                   momentum_stable = momentum_stable,
                   nesterov_accelerated_gradient = nesterov_accelerated_gradient,
                   input_dropout_ratio = input_dropout_ratio,
                   hidden_dropout_ratios = hidden_dropout_ratios,
                   l1 = l1,
                   l2 = l2,
                   max_w2 = max_w2,
                   initial_weight_distribution = match.arg(initial_weight_distribution),
                   initial_weight_scale = initial_weight_scale,
                   initial_weights = initial_weights,
                   initial_biases = initial_biases,
                   loss = match.arg(loss),
                   distribution = distribution,
                   quantile_alpha = quantile_alpha,
                   tweedie_power = tweedie_power,
                   huber_alpha = huber_alpha,
                   score_interval = score_interval,
                   score_training_samples = score_training_samples,
                   score_validation_samples = score_validation_samples,
                   score_duty_cycle = score_duty_cycle,
                   classification_stop = classification_stop,
                   regression_stop = regression_stop,
                   stopping_rounds = stopping_rounds,
                   stopping_metric = match.arg(stopping_metric),
                   stopping_tolerance = stopping_tolerance,
                   max_runtime_secs = max_runtime_secs,
                   score_validation_sampling = match.arg(score_validation_sampling),
                   diagnostics = diagnostics,
                   fast_mode = fast_mode,
                   force_load_balance = force_load_balance,
                   variable_importances = variable_importances,
                   replicate_training_data = replicate_training_data,
                   single_node_mode = single_node_mode,
                   shuffle_training_data = shuffle_training_data,
                   missing_values_handling = match.arg(missing_values_handling),
                   quiet_mode = quiet_mode,
                   autoencoder = autoencoder,
                   sparse = sparse,
                   col_major = col_major,
                   average_activation = average_activation,
                   sparsity_beta = sparsity_beta,
                   max_categorical_features = max_categorical_features,
                   reproducible = reproducible,
                   export_weights_and_biases = export_weights_and_biases,
                   mini_batch_size = mini_batch_size,
                   categorical_encoding = match.arg(categorical_encoding),
                   elastic_averaging = elastic_averaging,
                   elastic_averaging_moving_rate = elastic_averaging_moving_rate,
                   elastic_averaging_regularization = elastic_averaging_regularization) 
}
