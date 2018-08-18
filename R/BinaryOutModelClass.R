logit_linkinv <- function(x) { plogis(x) }

#----------------------------------------------------------------------------------
# Classes for predicting based on fitted regression models with binary outcome Bin ~ Xmat
#----------------------------------------------------------------------------------
# predict_single_reg <- function(self) UseMethod("predict_single_reg")
predict_single_reg <- function(self) {
  fitfunname <- self$getfit$fitfunname
  if (fitfunname %in% c("speedglm", "glm")) { return(predict_single_reg.glm(self)) }
  if (fitfunname == "h2o.ensemble") { return(predict_single_reg.h2o(self)) }
  if (fitfunname == "sl3_pipelines") { return(predict_single_reg.sl3(self)) }
  if (fitfunname == "SuperLearner") { return(predict_single_reg.SL(self)) }
}

predict_single_reg.long <- function(self) {
  fitfunname <- self$getfit$fitfunname
  if (fitfunname %in% c("speedglm", "glm")) { return(predict_single_reg.glm.long(self)) }
  if (fitfunname == "h2o.ensemble") { return(predict_single_reg.h2o(self)) }
  if (fitfunname == "SuperLearner") { return(predict_single_reg.SL(self)) }
}

# Generic prediction fun for logistic regression coefs, predicts P(A = 1 | newXmat)
predict_single_reg.glm <- function(self) {  # Does not handle cases with deterministic Anodes in the original data.
  model.fit <- self$getfit
  Xmat <- self$getXmat  
  assert_that(!is.null(Xmat)); assert_that(!is.null(self$subset_idx))
  pAout <- rep.int(gvars$misval, self$n)  # Set to default missing value for A[i] degenerate/degerministic/misval
  if (sum(self$subset_idx > 0)) {
    eta <- Xmat[,!is.na(model.fit$coef), drop = FALSE] %*% model.fit$coef[!is.na(model.fit$coef)]
    pAout[self$subset_idx] <- match.fun(FUN = model.fit$linkfun)(eta)
  }
  return(pAout)
}

predict_single_reg.glm.long <- function(self) {
  ID <- NULL  # To fix problem of no visible binding for global variable "ID"
  model.fit <- self$getfit
  Xmat <- self$getXmat  
  Y_vals <- self$getY
  assert_that(!is.null(Xmat)); assert_that(nrow(Xmat)==length(Y_vals)); assert_that(!is.null(self$subset_idx))
  pAout <- rep.int(gvars$misval, self$n)
  if (sum(self$subset_idx > 0)) {
    # -----------------------------------------------------------------
    # OBTAINING PREDICTIONS FOR LONG FORMAT P(Ind_j = 1 | Bin_j, W) BASED ON EXISTING POOLED FIT:
    # -----------------------------------------------------------------
    eta <- Xmat[,!is.na(model.fit$coef), drop = FALSE] %*% model.fit$coef[!is.na(model.fit$coef)]
    probA1 <- match.fun(FUN = model.fit$linkfun)(eta)
    # -----------------------------------------------------------------
    # GETTING ID-BASED PREDICTIONS (n) as cumprod of P(Ind_j = 1 | Bin_j, W) for j = 1, ..., K
    # -----------------------------------------------------------------    
    ProbAeqa_long <- as.vector(probA1^(Y_vals) * (1L - probA1)^(1L - Y_vals))    
    res_DT <- data.table(ID = self$ID, ProbAeqa_long = ProbAeqa_long)
    res_DT <- res_DT[, list(cumprob = cumprod(ProbAeqa_long)), by = ID]
    data.table::setkeyv(res_DT, c("ID")) # sort by ID
    res_DT_short <- res_DT[unique(res_DT[, key(res_DT), with = FALSE]), mult = 'last']
    ProbAeqa <- res_DT_short[["cumprob"]]
    pAout[self$subset_idx] <- ProbAeqa
  }
  return(pAout)
}

predict_single_reg.h2o <- function(self) {
  model.fit <- self$getfit$model.fit
  Xmat <- self$getXmat  
  assert_that(!is.null(Xmat)); assert_that(!is.null(self$subset_idx))
  pAout <- rep.int(gvars$misval, self$n)
  if ( any(class(model.fit) %in% "h2o.ensemble")) {
    if (sum(self$subset_idx > 0)) {
      test <- h2o::as.h2o(Xmat)
      predictions <- predict(model.fit, newdata = test)
      if (model.fit$family == "gaussian") { 
        pAout[self$subset_idx] <-  as.vector(predictions$pred)
      } else if (model.fit$family == "binomial") {
        pAout[self$subset_idx] <-  as.vector(predictions$pred[, 3])
      }
    }
  }
  return(pAout)
}

predict_single_reg.SL <- function(self) {
  model.fit <- self$getfit$model.fit
  Xmat <- self$getXmat  
  assert_that(!is.null(Xmat)); assert_that(!is.null(self$subset_idx))
  pAout <- rep.int(gvars$misval, self$n)
  if ( any(class(model.fit) %in% "SuperLearner")) {
    if (sum(self$subset_idx > 0)) {
      test <- data.frame(Xmat)
      predictions <- predict(model.fit, newdata = test, onlySL = TRUE)
      pAout[self$subset_idx] <-  as.vector(predictions$pred)
    }
  }
  return(pAout)
}

predict_single_reg.sl3(self) {
  model.fit <- self$getfit$model.fit
  Xmat <- self$getXmat  
  assert_that(!is.null(Xmat)); assert_that(!is.null(self$subset_idx))
  pAout <- rep.int(gvars$misval, self$n)
  if ( any(class(model.fit) %in% "sl3")) {
    if (sum(self$subset_idx > 0)) {
      test <- data.frame(Xmat)
      predictions <- model.fit$predict()
      pAout[self$subset_idx] <-  as.vector(predictions$pred)
    }
  }
  return(pAout)
}

#----------------------------------------------------------------------------------
# Classes for fitting regression models with binary outcome Bin ~ Xmat
#----------------------------------------------------------------------------------
# fit_single_reg <- function(self) UseMethod("fit_single_reg")
fit_single_reg <- function(self) {
  estimator <- self$estimator
  if (estimator == "glm__glm") { return(fit_single_reg.glmS3(self)) }
  if (estimator == "speedglm__glm") { return(fit_single_reg.speedglmS3(self)) }
  if (estimator == "h2o__ensemble") { return(fit_single_reg.h2oS3(self)) } 
  if (estimator == "sl3_pipelines") { return(fit_single_reg.sl3S3(self)) }
  if (estimator == "SuperLearner") { return(fit_single_reg.SLS3(self)) }
}
      
# S3 method for glm binomial family fit, takes BinaryOutModel object:
fit_single_reg.glmS3 <- function(self) {
  Xmat <- self$getXmat
  Y_vals <- self$getY
  wt_vals <- self$getWeight
  if (gvars$verbose) {
    print("calling glm.fit...")
    print("number of observations: " %+% nrow(Xmat))
  }
  
  if (nrow(Xmat) == 0L) {  # Xmat has 0 rows: return NA's and avoid throwing exception:
    model.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    ctrl <- glm.control(trace=FALSE, maxit=1000)
    # scale weights because there were rare problems where large weights caused convergence problems
    model.fit <- stats::glm.fit(x = Xmat, y = Y_vals, weights = as.vector(scale(wt_vals, center = FALSE)), 
                                family = quasibinomial(), control = ctrl)
  }
  fit <- list(coef = model.fit$coef, linkfun = "logit_linkinv", fitfunname = "glm")
  if (gvars$verbose) print(fit$coef)
  class(fit) <- c(class(fit), c("glmS3"))
  return(fit)
}
        
# S3 method for speedglm binomial family fit, takes BinaryOutModel object:
fit_single_reg.speedglmS3 <- function(self) {
  Xmat <- self$getXmat
  Y_vals <- self$getY
  wt_vals <- self$getWeight
  if (gvars$verbose) {
    print("calling speedglm.wfit...")    
    print("number of observations: " %+% nrow(Xmat))
  }
        
  if (nrow(Xmat) == 0L) { # Xmat has 0 rows: return NA's and avoid throwing exception
    model.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
  } else {
    model.fit <- try(
      # scale weights because there were rare problems where large weights caused convergence problems
      speedglm::speedglm.wfit(X = Xmat, y = Y_vals, weights = as.vector(scale(wt_vals, center = FALSE)), 
                              family = quasibinomial()), silent = TRUE)
    if (inherits(model.fit, "try-error")) {  # if failed, fall back on stats::glm
      message("speedglm::speedglm.wfit failed, falling back on stats:glm.fit; ", model.fit)
      return(fit_single_reg.glmS3(self))
    } 
  }
  if (gvars$verbose) print(model.fit$coef)
  fit <- list(coef = model.fit$coef, linkfun = "logit_linkinv", fitfunname = "speedglm")
  class(fit) <- c(class(fit), c("speedglmS3"))
  return(fit)
}        
      
# S3 method for h2o.ensemble binomial family fit, takes BinaryOutModel object:
fit_single_reg.h2oS3 <- function(self) {
  Xmat <- self$getXmat
  Y_vals <- self$getY
  wt_vals <- self$getWeight  # Cannot implement it since no weight argument in h2oEnsemble
  if (length(unique(wt_vals)) > 1) {
    message("Choose other estimators if want to implement obs.wt since h2o.ensemble does NOT have weight argument.")
  }
  h2olearner <- getopt("h2olearner")
  h2ometalearner <- getopt("h2ometalearner")
  CVfolds <- getopt("CVfolds")
  if (gvars$verbose) {
    print("calling h2o.ensemble...")
    print(length(h2olearner) %+% " machine learning algorithm(s): " %+% paste0(h2olearner, collapse = '; '))
    print("h2o metalearner: " %+% h2ometalearner)
    print("number of observations: " %+% nrow(Xmat))
  } 
  
  if (nrow(Xmat) == 0L) { # Xmat has 0 rows or : return NA's and avoid throwing exception
    model.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("No observations fall into the bin" %+% " So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n")
    }
  } else if (nrow(Xmat) != 0L & length(unique(Y_vals)) == 1) {  # Responses are constant
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("Responses are constant. So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n" %+% "falling back on speedglm::speedglm.wfit;\n")
    }
    return(fit_single_reg.speedglmS3(self))
  } else {
    localH2O <- h2o::h2o.init(nthreads = -1)
    # h2o.removeAll()
    if (!all(grepl("glm", h2olearner))) {
      message("================================================================")
      message("Caution: h2olearner include other algorithms except GLM (and GLMLNET)\n"
              %+% "such as randomForest, GBM and deeplearning. These aggressive\n" %+%
                "algorithms may provide a different result than only using glm\n")
      message("================================================================") 
    }
    train <- data.frame(Xmat[, colnames(Xmat)[colnames(Xmat) != "Intercept"]], Y_vals)
    names(train)[ncol(train)] <- self$outvar
    train <- h2o::as.h2o(train)
    if (all(Y_vals %in% 0:1)) {
      train[, self$outvar] <- h2o::as.factor(train[, self$outvar])
      h2oFamily <- "binomial"
    } else {
      h2oFamily <- "gaussian"
    }
    x <- setdiff(names(train), self$outvar) 
    model.fit <- try(h2oEnsemble::h2o.ensemble(
      x = x, y = self$outvar, training_frame = train, family = h2oFamily,  
      learner = h2olearner, metalearner = h2ometalearner, cvControl = list(V=CVfolds)))
    if (inherits(model.fit, "try-error")) {  # if failed, fall back on SuperLearner::SuperLearner
      message("h2oEnsemble::h2o.ensemble failed, falling back on SuperLearner::SuperLearner; \n", model.fit)
      return(fit_single_reg.SLS3(self)) 
    } 
    if (gvars$verbose) {
      if (h2ometalearner == "h2o.glm.wrapper") {
        print("the coefficients of the model " %+% self$outvar); print(h2o::h2o.coef(model.fit$metafit))
        print("Note, the first value is for the intercept")
      } else {
        print("there is no coefficient here since the metalearner is not GLM")
      }
    }       
  }
  fit <- list(model.fit = model.fit, coef = NULL, h2o.library = h2olearner, fitfunname = "speedglm")
  if (class(model.fit) == "h2o.ensemble" & h2ometalearner == "h2o.glm.wrapper") fit$coef <- h2o::h2o.coef(model.fit$metafit)
  if (class(model.fit) == "h2o.ensemble") {
    class(fit) <- c(class(fit), "h2o.ensemble")
    fit$fitfunname <- "h2o.ensemble"
    fit$family <- h2oFamily
  }
  return(fit)
}

# S3 method for sl3 binomial family fit, takes BinaryOutModel objects:
fit_single_reg.sl3S3 <-  function(self) {
  Xmat <- self$getXmat
  Y_vals <- self$getY
  wt_vals <- self$getWeight
  learner <- getopt("sl3_learner")
  metalearner <- getopt("sl3_metalearner")
  if (gvars$verbose) {
    print("calling sl3...")
    print(length(SL.library) %+% " machine learning algorithm(s): " %+% paste0(names(learner), collapse = '; '))
    print("number of observations: " %+% nrow(Xmat))
  }    
  
  if (nrow(Xmat) == 0L) {  # Xmat has 0 rows or Responses are constant: return NA's and avoid throwing exception
    model.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("No observations fall into the bin. So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n")
    }
  } else if (nrow(Xmat) != 0L & length(unique(Y_vals)) == 1) {  # Responses are constant
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("Responses are constant. So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n" %+% "falling back on speedglm::speedglm.wfit;\n")
    }
    return(fit_single_reg.speedglmS3(self))
  } else {
    if (all(Y_vals >= 0 & Y_vals <= 1)) {sl3Family <- "binomial" } else { sl3Family <- "gaussian" }
    if (sl3Family == "binomial") {
      cat("#######################################################################\n")
      cat("Currently we only accept binomial outcome when using sl3_pipelines " %+% self$outvar %+% "\n" %+% "falling back on speedglm::speedglm.wfit;\n")
    }
    return(fit_single_reg.speedglmS3(self))
    
    n <- length(Y_vals)
    X <- data.frame(Xmat[, colnames(Xmat)[colnames(Xmat) != "Intercept"]])
    data <- cbind(X, y = Y_vals, weights = wt_vals)
    Wnodes <- names(X)
    Anode <- "y"
    task <- sl3::sl3_Task$new(data, covariates = Wnodes, outcome = Anode, weights = "weights")
    
    # define Super Learner
    binom_sl <- sl3::make_learner(Lrnr_sl, learners, logit_metalearner)
    
    model.fit <- try(binom_sl$train(task))
    
    if (inherits(model.fit, "try-error")) { # if failed, fall back on SuperLearner::SuperLearner
      message(" failed, falling back on SuperLearner::SuperLearner; ", model.fit)
      return(fit_single_reg.SLS3(self))
    }
  }
  fit <- list(model.fit = model.fit, coef = NULL, SL.library = SL.library, fitfunname = "speedglm")
  if (any(class(model.fit) == "Lrnr_sl")) {
    fit$coefficients <- model.fit$coefficients
    class(fit) <- c(class(fit), "sl3")
    fit$fitfunname <- "sl3"
  }
  return(fit)
}

# S3 method for SuperLearner binomial family fit, takes BinaryOutModel objects:
fit_single_reg.SLS3 <- function(self) {
  Xmat <- self$getXmat
  Y_vals <- self$getY
  wt_vals <- self$getWeight
  SL.library <- getopt("SL.library")
  CVfolds <- getopt("CVfolds")
  if (gvars$verbose) {
    print("calling SuperLearner...")
    print(length(SL.library) %+% " machine learning algorithm(s): " %+% paste0(SL.library, collapse = '; '))
    print("number of observations: " %+% nrow(Xmat))
  }    
  
  if (nrow(Xmat) == 0L) {  # Xmat has 0 rows or Responses are constant: return NA's and avoid throwing exception
    model.fit <- list(coef = rep.int(NA_real_, ncol(Xmat)))
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("No observations fall into the bin. So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n")
    }
  } else if (nrow(Xmat) != 0L & length(unique(Y_vals)) == 1) {  # Responses are constant
    if (gvars$verbose) {
      cat("#######################################################################\n")
      cat("Responses are constant. So NO MODEL FITTED for this stratum " %+% self$outvar %+% "\n" %+% "falling back on speedglm::speedglm.wfit;\n")
    }
    return(fit_single_reg.speedglmS3(self))
  } else {
    if (all(Y_vals >= 0 & Y_vals <= 1)) { SLFamily <- binomial() } else { SLFamily <- gaussian() }
    n <- length(Y_vals)
    X <- data.frame(Xmat[, colnames(Xmat)[colnames(Xmat) != "Intercept"]])
    model.fit <- try(SuperLearner::SuperLearner( 
      Y = Y_vals, X = X, family = SLFamily, verbose = gvars$verbose, SL.library = SL.library, 
      obsWeights = as.vector(scale(wt_vals, center = FALSE)), cvControl = list(V=CVfolds), 
      control = list(saveFitLibrary=TRUE)))
    
    if (inherits(model.fit, "try-error")) { # if failed, fall back on speedglm::speedglm.wfit
      message("SuperLearner::SuperLearner failed, falling back on speedglm::speedglm.wfit; ", model.fit)
      return(fit_single_reg.speedglmS3(self))
    }
  }
  fit <- list(model.fit = model.fit, coef = NULL, SL.library = SL.library, fitfunname = "speedglm")
  if (class(model.fit) == "SuperLearner") {
    fit$coef <- model.fit$coef
    class(fit) <- c(class(fit), "SuperLearner")
    fit$fitfunname <- "SuperLearner"
  }
  return(fit)
}

# S3 methods for getting coefs from fitted BinaryOutModel class object
coef.BinaryOutModel <- function(binaryoutmodel) {
  assert_that(binaryoutmodel$is.fitted)
  binaryoutmodel$getfit$coef
}

summary.BinaryOutModel <- function(binaryoutmodel) {
  assert_that(binaryoutmodel$is.fitted)
  fit <- binaryoutmodel$getfit
  append(list(reg = binaryoutmodel$show()), fit)
}
  
# -------------------------------------- binirized.to.DTlong --------------------------------
# Purpose: Convert existing Bin matrix (Bin indicators) for continuous self$outvar into long format data.table with 3 columns:
# ID - row number; sVar_allB.j - bin indicators collapsed into one col; bin_ID - bin number identify for prev. columns
# automatically removed all missing (degenerate) bin indicators      
# -------------------------------------------------------------------------------------------
binirized.to.DTlong <- function(BinsDat_wide, binID_seq, ID, bin_names, pooled_bin_name, name.sVar) {
  # Convert Bin matrix into a data.table (without data.frame as intermediate), with new row ID column:
  DT_BinsDat_wide <- data.table::as.data.table(BinsDat_wide)[, c("ID") := ID]
  data.table::setcolorder(DT_BinsDat_wide, c("ID", names(DT_BinsDat_wide)[-ncol(DT_BinsDat_wide)]))
  # melt into long format:
  # Only the bin where the obs falls into receives 1, others receive 0
  sVar_melt_DT <- melt(DT_BinsDat_wide,
                      id.vars = "ID",
                      measure.vars = bin_names,
                      value.name = pooled_bin_name,  
                      variable.name = name.sVar,
                      variable.factor = FALSE,
                      na.rm = FALSE)
  nbin_rep <- rep(binID_seq, each = nrow(BinsDat_wide))
  # 1) Add bin_ID; 2) remove a column with Bin names; 3) remove all rows with NA value for outcome (degenerate bins)
  if (!is.data.table(sVar_melt_DT)) {
    class(sVar_melt_DT)
    stop("sVar_melt_DT is not a data.table")
  }
  # Remove column named as name.sVar; Remove rows who have NA under column named as pooled_bin_name
  sVar_melt_DT <- sVar_melt_DT[, c("bin_ID") := list(nbin_rep)][, (name.sVar) := NULL][!is.na(get(pooled_bin_name))]
  data.table::setkeyv(sVar_melt_DT, c("ID", "bin_ID"))  # sort by ID, bin_ID to prepare for merge with predictors (W)
  return(sVar_melt_DT)
}
     
# ------------------------------------------- join.Xmat -------------------------------------
# Purpose: Prepare predictors (W/X_mat) as data.table, adding row IDs for a join
# Join with sVar_melt_DT that is already in long format; Check that created IDs match exactly for both datasets
# -------------------------------------------------------------------------------------------
join.Xmat = function(X_mat, sVar_melt_DT, ID) {
  nIDs <- length(unique(sVar_melt_DT[["ID"]]))
  assert_that(nIDs == nrow(X_mat))
  X_mat_DT <- data.table::as.data.table(X_mat)[, c("ID") := ID]
  data.table::setkeyv(X_mat_DT, c("ID")) # sort by ID
  sVar_melt_DT <- sVar_melt_DT[X_mat_DT] # Merge long format (self$pooled_bin_name, binIDs) with predictors (W)
  return(sVar_melt_DT)
}      
        

## ---------------------------------------------------------------------
#' R6 class for modeling (fitting and predicting) for a single binary regression model P(B | PredVars)
#'
#'  \code{BinaryOutModel} can store and manage the (binarize/ discretized) design matrix Xmat and the outcome Bin for the binary regression 
#'  P(Bin|Xmat). It provides argument \code{self$estimator} to include different candidate estimators in the fitting and predicting library,  
#'  such as data-adaptive super learner algorithms and parametric logistic regression. When fitting one pooled regression across multiple 
#'  bins, it provides method to convert data from wide to long format when requested (to gain computational efficiency). 
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{bin_names} - Character vector of names of the bins.
#' \item{ID} - Integer vector of observation IDs used for pooling. \code{1:n}.
#' \item{pooled_bin_name} - Original name of the continuous covariate that was discretized into bins and then pooled.
#' \item{nbins} - Number of bins used for estimation of a continuous outvar, defined in ContinModel$new().
#' \item{estimator} - Character, one of "speedglm__glm" (default), "glm__glm", "h2o__ensemble", "SuperLearner".
#' \item{outvar} - Character, outcome name.
#' \item{predvars} - Character vector of predictor names.
#' \item{cont.sVar.flag} - Logical. If TRUE, indicate the original outcome variable is continuous.
# \item{fitclass} - (NOT IMPLEMENTED) Fit class based on estimator type.
#' \item{bw.j} - Bin width of a bin indicator obtained from the discretization of a continous covariate.
#' \item{is.fitted} - Logical. If TRUE, indicate the \code{BinaryOutModel} class object is fitted already.
#' \item{pool_cont} - Logical. If TRUE, perform pooling of bins.
#' \item{outvars_to_pool} - Character vector of outcome bin names for pooling.
#' \item{ReplMisVal0} - Logical. If TRUE, user-supplied gvars$misXreplace (Default to 0) will be used to replace all gvars$misval 
#'   among predictors. \code{ReplMisVal0} in \code{RegressionClass} will be used when instantiating an new object of \code{BinaryOutModel}. 
#' \item{n} - Number of rows in the input data.
#' \item{subset_expr} - Vector of length \code{n} that specifies a subset of data to be used in the fitting process. 
#'   Either logical, expression or indices.
#' \item{subset_idx} - Logical version of \code{subset_expr}. 
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg)}}{Use \code{reg} (a \code{\link{RegressionClass}} class object) to instantiate an new object of \code{BinaryOutModel} 
#'     for a single binary regression.}
#'   \item{\code{newdata(newdata, getoutvar = TRUE, ...)}}{Evaluate subset and perform correct subseting of data to construct 
#'     X_mat, Yvals & wt_vals.}
#'   \item{\code{define.subset.idx(data)}}{Create a logical vector which is converted from subset_expr}
#'   \item{\code{fit(overwrite = FALSE, data, predict = FALSE, savespace = TRUE, ...)}}{fit a binary regression. Note that \code{overwrite} is 
#'     Logical. If \code{FALSE} (Default), the previous fitted model cannot be overwritten by new fitting model. \code{savespace} is Logical. 
#'     If \code{TRUE} (Default), wipe out all internal data when doing many stacked regressions.}
#'   \item{\code{copy.fit(bin.out.model)}}{Take fitted BinaryOutModel object as an input and save the fit to itself.}
#'   \item{\code{predict(newdata, savespace = TRUE, ...)}}{Predict the response P(A = 1|W = w, E = e).}
#'   \item{\code{copy.predict(bin.out.model)}}{Tke BinaryOutModel object that contains the predictions for P(A=1|w,e) and save to itself}
#'   \item{\code{predictAeqa(newdata, bw.j.sA_diff, savespace = TRUE, wipeProb = TRUE)}}{Predict the response P(A = a|W = w, E = e) for 
#'     observed A, W, E. Note that wipeProb is logical argument for self$wipe.alldat. If FALSE, vectors of probA1 & probAeqa will be kept.}
#'   \item{\code{show()}}{Print regression formula, including outcome and predictor names.}
#' } 
#' @section Active Bindings:
#' \describe{
#'   \item{\code{wipe.alldat(wipeProb = TRUE)}}{...}
#'   \item{\code{getfit}}{...}
#'   \item{\code{getprobA1}}{...}
#'   \item{\code{getprobAeqa}}{...}
#'   \item{\code{emptydata}}{...}
#'   \item{\code{emptyY}}{...}
#'   \item{\code{emptyWeight}}{...}
#'   \item{\code{emptySubset_idx}}{...}
#'   \item{\code{getXmat}}{...}
#'   \item{\code{getY}}{...}
#'   \item{\code{getWeight}}{...}
#' }
#' @seealso \code{\link{DatKeepClass}}, \code{\link{RegressionClass}}, \code{\link{tmleCom_Options}}
#' @example tests/examples/1_BinaryOutModel_examples.R
#' @export
BinaryOutModel  <- R6Class(classname = "BinaryOutModel",
  cloneable = TRUE, # changing to TRUE to make it easy to clone input h_g0/h_gstar model fits
  portable = TRUE,
  class = TRUE,
  public = list(
    bin_names = NULL,
    ID = NULL,
    pooled_bin_name = NULL,
    nbins = integer(),
    estimator = character(), 
    outvar = character(),   # Outcome name(s)
    predvars = character(), # Names of predictor vars
    cont.sVar.flag = logical(),
    # fitclass = character(), # Fit class based on estimator type
    bw.j = numeric(),
    is.fitted = FALSE,   
    pool_cont = logical(),
    outvars_to_pool = character(),
    ReplMisVal0 = logical(),
    n = NA_integer_,        # Number of rows in the input data
    subset_expr = NULL,     # PASS THE LOGICAL EXPRESSIONS TO self$subset WHICH WILL BE EVALUTED IN THE ENVIRONMENT OF THE data
    subset_idx = NULL,      # Logical vector of length n (TRUE = include the obs)
 
    initialize = function(reg, ...) {
      assert_that(is.string(reg$outvar))
      assert_that(is.character(reg$predvars))
      self$outvar <- reg$outvar
      self$predvars <- reg$predvars
      self$pool_cont <- reg$pool_cont
      self$outvars_to_pool <- reg$outvars_to_pool
      self$ReplMisVal0 <- reg$ReplMisVal0
      self$nbins <- reg$nbins
      self$subset_expr <- reg$subset_vars
      if (is.null(reg$subset_vars)) {self$subset_expr <- TRUE}
      assert_that(is.logical(self$subset_expr) || is.call(self$subset_expr) || is.character(self$subset_expr))
      self$estimator <- reg$estimator
      # self$fitclass <- getopt("fitclass")      
      
      # Get the bin width (interval length) for the current bin name self$outvar (for discretized continuous A only):
      self$cont.sVar.flag <- self$outvar %in% names(reg$intrvls.width)
      if (self$cont.sVar.flag) {
        intrvl.idx <- which(names(reg$intrvls.width) %in% self$outvar)
        if (length(intrvl.idx) > 1) stop("non-unique names for intrvls.width in RegressionClass")
        self$bw.j <- reg$intrvls.width[intrvl.idx]
      } else {
        self$bw.j <- 1L
      }
      invisible(self)
    },
    
    # ----------------------------------------- newdata ------------------------------------------- 
    # Purpose: Sets X_mat, Yvals, wt_vals, evaluates subset and performs correct subseting of data  
    # --------------------------------------------------------------------------------------------- 
    newdata = function(newdata, getoutvar = TRUE, ...) {
      assert_that(is.DatKeepClass(newdata))
      # CALL self$setdata.long() when: 1) self$pool_cont is TRUE & 2) more than one outvars_to_pool
      if (self$pool_cont && length(self$outvars_to_pool)>1) {
        private$setdata.long(data = newdata, ...)
      } else {
        private$setdata(data = newdata, getoutvar, ...)
      }
      invisible(self)
    },
  
    # ------------------------------------ define.subset.idx -------------------------------------- 
    # Purpose: return a logical vector which is converted from Subset subset_expr
    # --------------------------------------------------------------------------------------------- 
    define.subset.idx = function(data) {
      if (is.logical(self$subset_expr)) {
        subset_idx <- self$subset_expr
      } else if (is.call(self$subset_expr)) {
        subset_idx <- data$evalsubst(subset_exprs = self$subset_expr)
      } else if (is.character(self$subset_expr)) {
        subset_idx <- data$evalsubst(subset_vars = self$subset_expr)
      }
      assert_that(is.logical(subset_idx))
      if ((length(subset_idx) < self$n) && (length(subset_idx) > 1L)) {
        if (gvars$verbose) message("subset_idx has smaller length than self$n; repeating subset_idx p times, for p: " %+% data$p)
        subset_idx <- rep.int(subset_idx, data$p)
        if (length(subset_idx) != self$n) stop("BinaryOutModel$define.subset.idx: self$n is not equal to nobs*p!")
      }
      assert_that((length(subset_idx) == self$n) || (length(subset_idx) == 1L))
      return(subset_idx)
    },   

    fit = function(overwrite = FALSE, data, predict = FALSE, savespace = TRUE, ...) { # Move overwrite to a field? ... self$overwrite
      if (gvars$verbose) print("fitting the model: " %+% self$show())
      if (!overwrite) assert_that(!self$is.fitted) # do not allow overwrite of previous fitted model unless explicitely asked  
      self$newdata(newdata = data, ...) # populate bindat with X_mat & Y_vals & wt_vals
      private$model.fit <- fit_single_reg(self = self) 
      self$is.fitted <- TRUE
      # **********************************************************************
      # to save RAM space when doing many stacked regressions no longer predicting in fit:
      # **********************************************************************
      if (savespace) self$wipe.alldat
      invisible(self)
    },

    # take fitted BinaryOutModel object as an input and save the fits to itself
    copy.fit = function(bin.out.model) {
      assert_that("BinaryOutModel" %in% class(bin.out.model))
      private$model.fit <- bin.out.model$getfit
      self$is.fitted <- TRUE
      invisible(self)
    },

    # Predict the response P(Bin = 1|W = w);
    # Does not need to know the actual values of the binary outcome Bin to do prediction ($newdata(, getouvar = FALSE))
    # P(A[i]=1|W=w): uses private$model.fit to generate predictions for newdata:
    predict = function(newdata, savespace = TRUE, ...) {
      assert_that(self$is.fitted)
      if (missing(newdata)) { stop("must provide newdata for BinaryOutModel$predict()") }
      # re-populate self with new X_mat:
      self$newdata(newdata = newdata, getoutvar = FALSE, ...) 
      if (self$pool_cont && length(self$outvars_to_pool) > 1) {
        stop("BinaryOutModel$predict is not applicable to pooled regression, call BinaryOutModel$predictAeqa instead")
      } else {
        private$probA1 <- predict_single_reg(self = self)
      }
      if (savespace) self$emptydata  # Xmat is no longer needed, but subset, outvar & probA1 may be needed for private$probA1
      invisible(self)
    },

    # take BinaryOutModel object that contains the predictions for P(A=1|W) and save these predictions to self$
    copy.predict = function(bin.out.model) {
      assert_that("BinaryOutModel" %in% class(bin.out.model))
      assert_that(self$is.fitted)
      private$probA1 <- bin.out.model$getprobA1
    },

    # Predict the response P(Bin = b|W = w), which is returned invisibly;
    # Needs to know the values of b for prediction ($newdata(, getouvar = TRUE))
    # WARNING: This method cannot be chained together with methods that follow (s.a, class$predictAeqa()$fun())
    predictAeqa = function(newdata, bw.j.sA_diff, savespace = TRUE, wipeProb = TRUE) { 
      # P(A[i]=a|W=w) - calculating the likelihood for indA[i] (n vector of a's)
      if (missing(newdata) && !is.null(private$probAeqa)) { return(private$probAeqa) }
      assert_that(self$is.fitted)
      assert_that(!missing(newdata))
      self$newdata(newdata = newdata, getoutvar = TRUE) # populate bindat with new design matrix covars X_mat
      assert_that(is.logical(self$subset_idx))
      n <- newdata$nobs
      # obtain predictions (likelihood) for response on fitted data (from long pooled regression):
      if (self$pool_cont && length(self$outvars_to_pool) > 1) {
        probAeqa <- predict_single_reg.long(self = self) # overwrite probA1 with new predictions:
      } else {
        # get predictions for P(A[j]=1|W=newdata) from newdata:
        probA1 <- predict_single_reg(self = self)
        indA <- newdata$get.outvar(self$subset_idx, var = self$outvar) # Always a vector of 0/1
        assert_that(is.integerish(indA)) # check that obsdat.A is always a vector of of integers
        probAeqa <- rep.int(1L, n) # for missing, the likelihood is always set to P(A = a) = 1.
        assert_that(!any(is.na(probA1[self$subset_idx]))) # check that predictions P(A=1 | dmat) exist for all obs.
        probA1 <- probA1[self$subset_idx]
        # discrete version for the joint density:
        probAeqa[self$subset_idx] <- probA1^(indA) * (1 - probA1)^(1L - indA)
        # continuous version for the joint density:
        # probAeqa[self$subset_idx] <- (probA1^indA) * exp(-probA1)^(1 - indA)
        # Alternative intergrating the last hazard chunk up to x:
        # difference of A value and its left most bin cutoff: x - b_{j-1}
        if (!missing(bw.j.sA_diff)) {
          # + integrating the constant hazard all the way up to value of each sa:
          # probAeqa[self$subset_idx] <- probAeqa[self$subset_idx] * (1 - bw.j.sA_diff[self$subset_idx]*(1/self$bw.j)*probA1)^(indA)
          # cont. version of above:
          probAeqa[self$subset_idx] <- probAeqa[self$subset_idx] * exp(-bw.j.sA_diff[self$subset_idx]*(1/self$bw.j)*probA1)^(indA)
        }
        private$probA1 <- probA1
      }
      private$probAeqa <- probAeqa
      # **********************************************************************
      # to save RAM space when doing many stacked regressions wipe out all internal data:
      if (savespace) { self$wipe.alldat <- wipeProb } #  if wipeProb = FALSE, probA1 & probAeqa will be kept
      # **********************************************************************
      return(probAeqa)
    },
      
    show = function() { "P(" %+% self$outvar %+% "|" %+% paste(self$predvars, collapse=",") %+% ")" }
  ),
        
  active = list(
    wipe.alldat = function(wipeProb = TRUE) {
      if (wipeProb) { private$probA1 <- NULL }
      if (wipeProb) { private$probAeqa <- NULL }
      self$emptydata
      self$emptyY
      # self$emptyWeight
      self$emptySubset_idx
      return(self)
    },
    getfit = function() { private$model.fit },
    getprobA1 = function() { private$probA1 },
    getprobAeqa = function() { private$probAeqa },    
    emptydata = function() { private$X_mat <- NULL },
    emptyY = function() { private$Y_vals <- NULL},
    emptySubset_idx = function() { self$subset_idx <- NULL },
    emptyWeight = function() { private$wt_vals <- NULL },
    getXmat = function() { private$X_mat },
    getY = function() { private$Y_vals },
    getWeight = function() { private$wt_vals }
  ),
        
  private = list(
    model.fit = list(),   # the model fit (either coefficients or the model fit object)
    probA1 = NULL,        # Predicted probA=1 conditional on X_mat
    probAeqa = NULL,      # Likelihood of observing a particular value A=a conditional on X_mat
    X_mat = NULL,
    Y_vals = NULL,
    wt_vals = NULL,
    
    setdata = function(data, getoutvar, ...) {
      # everything is performed using data$ methods (data is of class DatKeepClass)
      assert_that(is.DatKeepClass(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset.idx(data)
      private$wt_vals <- data$get.obsweights(rowsubset = self$subset_idx)
      if (getoutvar) private$Y_vals <- data$get.outvar(rowsubset = self$subset_idx, var = self$outvar) # Always a vector
      if (sum(self$subset_idx) == 0L) {  # When nrow(X_mat) == 0L avoids exception (when nrow == 0L => prob(A=a) = 1)
        private$X_mat <- matrix(, nrow = 0L, ncol = (length(self$predvars) + 1))
        colnames(private$X_mat) <- c("Intercept", self$predvars)
      } else {
        # *** THIS IS THE ONLY LOCATION IN THE PACKAGE WHERE CALL TO DatKeepClass$get.dat.sVar() IS MADE ***
        private$X_mat <- as.matrix(cbind(Intercept = 1, data$get.dat.sVar(self$subset_idx, self$predvars)))        
        if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace  # To find and replace misvals in X_mat
      }
      invisible(self)
    },
        
    setdata.long = function(data, ...) {
      assert_that(is.DatKeepClass(data))
      self$n <- data$nobs
      self$subset_idx <- self$define.subset.idx(data)
      if (!(data$active.bin.sVar %in% self$outvar)) { stop("currently binirized sVar does not match self$outvar argument") }
      # Setting up object fields related to pooling of continuous A:
      self$pooled_bin_name <- data$pooled.bin.nm.sVar(self$outvar)
      self$bin_names <- self$outvars_to_pool

      if (gvars$verbose) {
        print("self$bin_names: "); print(self$bin_names)
        print("self$pooled_bin_name: "); print(self$pooled_bin_name)
        print("self$outvar: "); print(self$outvar)
        print("self$nbins: "); print(self$nbins)
      }

      binID_seq <- 1L : self$nbins
      BinsDat_wide <- data$get.dat.sVar(self$subset_idx, covars = self$outvars_to_pool)
      self$ID <- as.integer(1:nrow(BinsDat_wide))

      # To grab bin Ind mat directly (prob a bit faster): BinsDat_wide <- data$dat.bin.sVar[self$subset_idx, ]
      BinsDat_long <- binirized.to.DTlong(BinsDat_wide = BinsDat_wide, binID_seq = binID_seq, ID = self$ID,
                                          bin_names = self$bin_names, pooled_bin_name = self$pooled_bin_name,
                                          name.sVar = self$outvar)
      X_mat.addwt <- cbind(data$get.dat.sVar(self$subset_idx, covars = self$predvars), weights = data$get.obsweights(self$subset_idx))
      sVar_melt_DT <- join.Xmat(X_mat = X_mat.addwt, sVar_melt_DT = BinsDat_long, ID = self$ID)
      # prepare design matrix for modeling w/ glm.fit or speedglm.wfit:
      # select bin_ID + predictors, add intercept column
      X_mat <- sVar_melt_DT[, c("bin_ID", self$predvars), with=FALSE][, c("Intercept") := 1] 
      setcolorder(X_mat, c("Intercept", "bin_ID", self$predvars)) # re-order columns by reference (no copy)
      self$ID <- sVar_melt_DT[["ID"]]
      private$X_mat <- as.matrix(X_mat)
      private$Y_vals <- sVar_melt_DT[, self$pooled_bin_name, with = FALSE][[1]] # outcome vector
      private$wt_vals <- sVar_melt_DT[, "weights", with = FALSE][[1]] # observation weight vector
      if (self$ReplMisVal0) private$X_mat[gvars$misfun(private$X_mat)] <- gvars$misXreplace
      if (gvars$verbose) {
        print("private$X_mat[1:10,]"); print(private$X_mat[1:10,])
        print("head(private$Y_vals)"); print(head(private$Y_vals, 100))
        print("head(private$wt_vals)"); print(head(private$wt_vals, 100))
      }
      invisible(self)
    } 
  )
)
