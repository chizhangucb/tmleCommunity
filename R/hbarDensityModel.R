#-----------------------------------------------------------------------------
# Methods for fitting and predicting the IPTW (clever covariate) for covariates (A, W, E): P_{g^*}(A | W, E)/P_{g0}(A | W, E)
# 2 ancillary functions: predict.hbars, fit.hbars
#-----------------------------------------------------------------------------
fitGenericDensity <- function(data, Anodes, Wnodes, gform = NULL, f_gstar, h.gstar_GenericModel = NULL, 
                              lbound = 0.025, n_MCsims = 1, obs.wts = NULL, 
                              rndseed = NULL, verbose = TRUE) {
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  gvars$verbose <- verbose
  message("Running fitGenericDensity with the following settings from tmleCom_Options(): "); str(gvars$opts)
  
  if (is.null(obs.wts)) obs.wts <- rep(1, nrow(data))
  
  ## Check if any unexpected inputs
  nodes <- list(Anodes = Anodes, Wnodes = Wnodes, Crossnodes = NULL)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, NULL, gform, NULL, "logistic", NULL, obs.wts)) stop()
  
  #----------------------------------------------------------------------------------
  # DEFINING (OPTIONAL) REGRESSION FORMS 
  #----------------------------------------------------------------------------------
  h.gstar.sVars <- define_regform(as.formula(gform), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes$Wnodes)
  
  if (verbose) {
    print("Input regression gform (P(A|W) under g.star): " %+% gform)
    print("Derived regression gform (P(A|W) under g.star): "); str(h.gstar.sVars)
  }
  
  ## Create data based on gform, in case of interaction or higher-order term.
  if (!is.null(gform)) {
    allcovRHS <- strsplit(deparse(as.formula(gform)[[3]]), " \\+ ")[[1]]
    merged.form <- reformulate(allcovRHS, response = NULL)  # Reformulate a formula including all legitimate character of the RHS 
    if (any(!allcovRHS %in% unique(c(unlist(nodes), names(data))))) {
      ExtraDat <- as.data.frame(model.matrix(merged.form, data = data))
      data <- cbind(data, ExtraDat[, setdiff(names(ExtraDat), c(names(data), "(Intercept)")), drop = FALSE])
      Extradat <- NULL
    }
  } else {
    merged.form <- NULL
  }
  nodes <- append(nodes, list(Crossnodes = setdiff(names(data), Reduce(c, nodes))))
  
  ## Create an R6 object that stores and manages the input data, later passed on to estimation algorithm(s)
  OData.ObsP0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
  OData.ObsP0$addObsWeights(obs.wts = obs.wts)
  
  #----------------------------------------------------------------------------------
  # Defining and estimating treatment mechanism P(A|W)
  #----------------------------------------------------------------------------------
  Ws.gstar_nms <- h.gstar.sVars$predvars
  As_nms_gstar <- h.gstar.sVars$outvars
  subsets_expr <- lapply(As_nms_gstar, function(var) {var})  # subsetting by !gvars$misval on A
  As_classes <- OData.ObsP0$type.sVar[As_nms_gstar]
  
  if (verbose) {
    message("================================================================")
    message("fitting h_gstar with covariates: ", "P(" %+% paste(As_nms_gstar, collapse = ",") 
            %+% " | " %+% paste(Ws.gstar_nms, collapse = ",") %+% ")")
    message("================================================================")
  }
  
  OData.gstar <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
  OData.gstar$addObsWeights(obs.wts = obs.wts)
  OData.gstar$make.dat.sVar(p = n_MCsims, f.g_fun = f_gstar, regform = merged.form)
  
  if (gvars$verbose) {
    print("Generated new exposure covariates by sampling A from f_gstar (OData.gstar): "); print(class(OData.gstar$dat.sVar))
    print(dim(OData.gstar$dat.sVar)); print(head(OData.gstar$dat.sVar));
  }
  
  regclass.gstar <- RegressionClass$new(outvar.class = As_classes,
                                        outvar = As_nms_gstar,
                                        predvars = Ws.gstar_nms,
                                        subset_vars = subsets_expr,
                                        estimator = getopt("gestimator"))
  # Define Intervals Under g_star to Be The Same as under g0:  
  genericmodels.gstar <- GenericModel$new(reg = regclass.gstar, DatKeepClass.g0 = OData.ObsP0)
  
  if (!is.null(h.gstar_GenericModel)) {
    # 1) deep copy the object with model fits to genericmodels.gstar
    genericmodels.gstar <- h.gstar_GenericModel$clone(deep = TRUE)
  } else {
    genericmodels.gstar$fit(data = OData.gstar)
  }
  
  h_gstar <- genericmodels.gstar$predictAeqa(newdata = OData.ObsP0)
  h_gstar[is.nan(h_gstar)] <- 0     # 0/0 detection
  h_gstar <- bound(h_gstar, c(0, 1/lbound))
  
  model.h.fit <- list(genericmodels.gstar = genericmodels.gstar, lbound = lbound)
  return(list(h_gstar = h_gstar, model.h.fit = model.h.fit, OData.gstar = OData.gstar))
}


# @title Predict h weights under g_0 and g_star using existing model.h.fit model fit
# @name pred.hbars
# @export
# fit models for m_gAi
predict.hbars <- function(newdatnet = NULL, model.h.fit) {
  lbound <- model.h.fit$lbound
  # use original fitted data for prediction
  if (is.null(newdatnet)) {
    stop("newdatnet argument must be not null; this feature is not implemented")
  }
  # PASS ENTIRE newdatnet which will get subset, rather than constructing cY_mtx...
  h_gN <- model.h.fit$genericmodels.g0$predictAeqa(newdata = newdatnet)
  h_gstar <- model.h.fit$genericmodels.gstar$predictAeqa(newdata = newdatnet)
  h_gstar_h_gN <- h_gstar / h_gN
  h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
  h_gstar_h_gN <- bound(h_gstar_h_gN, c(0,1/lbound))
  return(h_gstar_h_gN)
}

# @title Defining and fitting the clever covariate h under g_0 and g_star, i.e. models P(A[j] | W, E, A[j-1])
# @name fit.hbars
# @importFrom assertthat assert_that is.count
# @export
# fit models for m_gAi
fit.hbars <- function(OData.ObsP0, est_params_list) {
  data <- est_params_list$data
  nodes <- est_params_list$nodes
  obs.wts <- est_params_list$obs.wts
  lbound <- est_params_list$lbound
  merged.form <- est_params_list$merged.form
  estnames <- est_params_list$estnames
  n_MCsims <- est_params_list$n_MCsims
  h.g0.sVars <- est_params_list$h.g0.sVars
  h.gstar.sVars <- est_params_list$h.gstar.sVars
  h.g0_GenericModel <- est_params_list$h.g0_GenericModel
  h.gstar_GenericModel <- est_params_list$h.gstar_GenericModel
  f.g0 <- est_params_list$f.g0
  f.gstar <- est_params_list$f.gstar
  savetime <- est_params_list$savetime.fit.hbars
  
  if (!is.null(h.g0_GenericModel)) {
    message("NOTE: Predictions for P(A|W,E) under g0 will be based on the fitted model in h.g0_GenericModel, " %+%
            "and all modeling settings will be ignored")
    # verify h.g0_GenericModel is consistent with GenericModel 
    assert_that(inherits(h.g0_GenericModel, "GenericModel"))
  }
  h_gstar_SummariesModel <- est_params_list$h_gstar_SummariesModel
  if (!is.null(h.gstar_GenericModel)) {
    message("NOTE: Predictions for P(A^*|W,E) under f_gstar will be based on the fitted model in h.gstar_GenericModel, " %+%
            "and all modeling settings will be ignored")
    # verify h.gstar_GenericModel is consistent with GenericModel 
    assert_that(inherits(h.gstar_GenericModel, "GenericModel"))
  }
   
  if (!is.null(f.gstar)  | !savetime | ("TMLE_A" %in% estnames)) {
    #---------------------------------------------------------------------------------
    # Getting OBSERVED W & A 
    #---------------------------------------------------------------------------------
    # Baseline variable names / expressions
    Ws.g0_nms <- h.g0.sVars$predvars
    Ws.gstar_nms <- h.gstar.sVars$predvars
    
    # Exposure variable names / expressions:
    As_nms_g0 <- h.g0.sVars$outvars
    As_nms_gstar <- h.gstar.sVars$outvars
    if (!all(As_nms_g0 == As_nms_gstar)) 
      stop("the outcome variable names defined by regressions hform.g0 & hform.gstar are not identical;" 
           %+% " current implementation requires these to be the same.")
    subsets_expr <- lapply(As_nms_g0, function(var) {var})  # subsetting by !gvars$misval on A

    # As' class parameters:
    As_classes <- OData.ObsP0$type.sVar[As_nms_g0]
    
    
#-------- Stage 1 --------- Calculate h_gN  
    if (gvars$verbose) {
      message("================================================================")
      message("fitting h_g0 with covariates: ", "P(" %+% paste(As_nms_g0, collapse = ",") 
              %+% " | " %+% paste(Ws.g0_nms, collapse = ",") %+% ")")
      message("================================================================")
    }
    
    p_h0 <- ifelse(is.null(f.g0), 1, n_MCsims)
    # If !is.null(f.g_fun) then OData.g0$dat.sVar IS NOT THE OBSERVED data (sVar), but rather sVar data sampled under known g0.
    if (!is.null(f.g0)) {
      if (gvars$verbose) message("generating OData.g0 under known g0")
      OData.g0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
      OData.g0$addObsWeights(obs.wts = obs.wts)
      OData.g0$make.dat.sVar(p = p_h0, f.g_fun = f.g0, regform = merged.form)
      print("head(OData.g0$dat.sVar): "); print(head(OData.g0$dat.sVar))
    } else {
      OData.g0 <- OData.ObsP0
    }
    
    regclass.g0 <- RegressionClass$new(outvar.class = As_classes,
                                       outvar = As_nms_g0,
                                       predvars = Ws.g0_nms,
                                       subset_vars = subsets_expr,
                                       estimator = getopt("gestimator"))
    genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
    if (!is.null(h.g0_GenericModel)) {
      # 1) deep copy model fits in h_g0_SummariesModel to genericmodels.g0
      genericmodels.g0 <- h.g0_GenericModel$clone(deep = TRUE)
    } else {
      genericmodels.g0$fit(data = OData.g0)
    }
    h_gN <- genericmodels.g0$predictAeqa(newdata = OData.ObsP0)
    
#-------- Stage 2 --------- Calculate h_gstar  
    if (gvars$verbose) {
      message("================================================================")
      message("fitting h_gstar with covariates: ", "P(" %+% paste(As_nms_gstar, collapse = ",") 
              %+% " | " %+% paste(Ws.gstar_nms, collapse = ",") %+% ")")
      message("================================================================")
    }
    
    OData.gstar <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
    OData.gstar$addObsWeights(obs.wts = obs.wts)
    OData.gstar$make.dat.sVar(p = n_MCsims, f.g_fun = f.gstar, regform = merged.form)
    
    if (gvars$verbose) {
      print("Generated new exposure covariates by sampling A from f_gstar (OData.gstar): "); print(class(OData.gstar$dat.sVar))
      print(dim(OData.gstar$dat.sVar)); print(head(OData.gstar$dat.sVar));
    }
    
    regclass.gstar <- RegressionClass$new(outvar.class = As_classes,
                                          outvar = As_nms_gstar,
                                          predvars = Ws.gstar_nms,
                                          subset_vars = subsets_expr,
                                          estimator = getopt("gestimator"))
    # Define Intervals Under g_star to Be The Same as under g0:  
    genericmodels.gstar <- GenericModel$new(reg = regclass.gstar, DatKeepClass.g0 = OData.g0)
    # Define Intervals Under g_star Based on exposure covariates Generated under g_star:
    # genericmodels.gstar <- GenericModel$new(reg = regclass.gstar, DatKeepClass.g0 = OData.gstar)
    
    if (!is.null(h.gstar_GenericModel)) {
      # 1) deep copy the object with model fits to genericmodels.gstar
      genericmodels.gstar <- h.gstar_GenericModel$clone(deep = TRUE)
    } else {
      genericmodels.gstar$fit(data = OData.gstar)
    }
    h_gstar <- genericmodels.gstar$predictAeqa(newdata = OData.ObsP0)
    
#-------- Stage 3 ---------  Calculate final h_bar (h_tilde) as ratio of h_gstar / h_gN and bound it
    h_gstar_h_gN <- h_gstar / h_gN
    h_gstar_h_gN[is.nan(h_gstar_h_gN)] <- 0     # 0/0 detection
    h_gstar_h_gN <- bound(h_gstar_h_gN, c(0, 1/lbound))
  
  } else {  # is.null(f.gstar) & savetime = TRUE & !("TMLE_A" %in% estnames)
    if (gvars$verbose)  {
      message("================================================================")
       message("Since f.gstar = NULL & using tmle.intercept in targeting step, in order to save time, \n"
               %+% "we skip g0 & gstar fitting procedure and directly set h_gstar_h_gN = 1 for each observation")
      message("================================================================")
    }
    
    OData.gstar <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
    OData.gstar$make.dat.sVar(p = n_MCsims, f.g_fun = f.gstar, regform = merged.form)
    if (gvars$verbose) {
      print("Generated new exposure covariates by sampling A from f_gstar (OData.gstar): "); print(class(OData.gstar$dat.sVar))
      print(dim(OData.gstar$dat.sVar)); print(head(OData.gstar$dat.sVar));
    }    
    
    if (!is.null(h.g0_GenericModel)) {
      genericmodels.g0 <- h.g0_GenericModel$clone(deep = TRUE)
    } else {
      genericmodels.g0 <- NULL
    }
    
    if (!is.null(h.gstar_GenericModel)) {      
      genericmodels.gstar <- h.gstar_GenericModel$clone(deep = TRUE)
    } else {
      genericmodels.gstar <- NULL
    }
   
    h_gstar_h_gN <- rep_len(1, OData.gstar$nobs)
  }
  model.h.fit <- list(genericmodels.g0 = genericmodels.g0, genericmodels.gstar = genericmodels.gstar, lbound = lbound)
  return(list(h_gstar_h_gN = h_gstar_h_gN, model.h.fit = model.h.fit, OData.gstar = OData.gstar))
}
