tmleSingleCommunity <- function(data, Ynode, Anodes, Wnodes, Enodes = NULL, YnodeDet = NULL, 
                                f_gstar1, f_gstar2 = NULL, Qform = NULL, Qbounds = NULL, alpha = 0.995, fluctuation = "logistic",                                                     
                                f.g0 = NULL, hform.g0 = NULL, hform.gstar = NULL, lbound = 0.005, obs.wts = NULL, 
                                h.g0_GenericModel = NULL, h.gstar_GenericModel = NULL, savetime.fit.hbars = TRUE, 
                                TMLE.targetStep = c("tmle.intercept", "tmle.covariate"),
                                n_MCsims = 1, 
                                CI_alpha = 0.05, 
                                rndseed = NULL, 
                                verbose = TRUE) {
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  gvars$verbose <- verbose
  message("Running tmleCommunity with the following settings from tmleCom_Options(): "); str(gvars$opts)
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (is.null(savetime.fit.hbars)) savetime.fit.hbars <- getopt("savetime.fit.hbars")
  if (is.null(obs.wts)) obs.wts <- rep(1, nrow(data)) 
  if (!is.null(Qform) && !is.null(Ynode)) {
    Qform <- paste(Ynode, substring(Qform, first = as.numeric(gregexpr("~", Qform))))
    message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+% Ynode)
    message("Thus the Qform becomes " %+% Qform)
  }
  if (!is.null(Qform) && is.null(Ynode)) {
    Ynode <- LhsVars(Qform)[1]
    message("Setting the Ynode to: " %+% Ynode)
  }
  TMLE.targetStep <- TMLE.targetStep[1]
  
  ## Check if any unexpected inputs
  if (!(TMLE.targetStep %in% c("tmle.intercept", "tmle.covariate"))) 
    stop("TMLE.targetStep argument must be either 'tmle.intercept' or 'tmle.covariate'")
  nodes <- list(Ynode = Ynode, Anodes = Anodes, Wnodes = Wnodes, Enodes = Enodes, Crossnodes = NULL)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, Qform, hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts)) stop()
  maptoYstar <- fluctuation=="logistic"  # if TRUE, cont Y values shifted & scaled to fall b/t (0,1)
  
  #----------------------------------------------------------------------------------
  # DEFINING (OPTIONAL) REGRESSION FORMS 
  #----------------------------------------------------------------------------------
  Q.sVars <- define_regform(as.formula(Qform), Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "Wnodes", "Enodes")])
  h.g0.sVars <- define_regform(as.formula(hform.g0), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[3:5])
  if (!is.null(hform.gstar)) {
    h.gstar.sVars <- define_regform(as.formula(hform.gstar), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[c("Wnodes", "Enodes")])
  } else {
    h.gstar.sVars <- h.g0.sVars
  }
  
  if (verbose) {
    print("Input regression Qform (E(Y|A,W,E)): " %+% Qform)
    print("Derived regression Qform (E(Y|A,W,E)):"); str(Q.sVars)
    print("Input regression hform.g0 (P(A|W,E) under g0): " %+% hform.g0)
    print("Derived regression hform.g0 (P(A|W,E) under g0): "); str(h.g0.sVars)
    print("Input regression hform.gstar (P(A|W,E) under g.star): " %+% hform.gstar)
    print("Derived regression hform.gstar (P(A|W,E) under g.star): "); str(h.gstar.sVars)
  }
  
  ## Create data based on Qform, hform.g0 and hform.gstar, in case of interaction or higher-order term.
  if (!is.null(c(Qform, hform.g0, hform.gstar))) {
    allcovRHS <- unique(unlist(lapply(c(Qform, hform.g0, hform.gstar), FUN = function(x) { strsplit(deparse(as.formula(x)[[3]]), " \\+ ")[[1]] })))
    merged.form <- reformulate(allcovRHS, response = NULL)  # Reformulate a formula including all legitimate character of the RHS in 3 formulae
    if (any(!allcovRHS %in% unique(c(unlist(nodes), names(data))))) {
      ExtraDat <- as.data.frame(model.matrix(merged.form, data = data))
      data <- cbind(data, ExtraDat[, setdiff(names(ExtraDat), c(names(data), "(Intercept)")), drop = FALSE])
      ExtraDat <- NULL
    }
  } else {
    merged.form <- NULL
  }
  nodes <- append(nodes, list(Crossnodes = setdiff(names(data), Reduce(c, nodes))))
  
  ## Create an R6 object that stores and manages the input data, later passed on to estimation algorithm(s)
  inputYs <- CreateInputs(data[, Ynode], Qbounds, alpha, maptoYstar)
  data[, Ynode] <- inputYs$Ystar
  OData.ObsP0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
  OData.ObsP0$addYnode(YnodeVals = inputYs$Ystar)
  OData.ObsP0$addObsWeights(obs.wts = obs.wts)
  nobs <- OData.ObsP0$nobs
  if (is.null(YnodeDet)) {
    determ.Q <- rep_len(FALSE, nobs)
  } else {
    determ.Q <- (data[, YnodeDet] == 1)
  }
  if (length(unique(obs.wts)) > 1 && any(unlist(OData$type.sVar[Anodes]) != "binary")) {
    warning("obs.wts are currently implemented on binary A. The results for non-binary A with weights may be unrealiable.")
  }
  
  #----------------------------------------------------------------------------------
  # Defining and estimating outcome mechanism E(Y|A, E, W)
  #----------------------------------------------------------------------------------
  if (verbose) {
    message("================================================================")
    message("fitting E(Y|A,W,E):= ", "P(" %+% nodes$Ynode %+% "=1 | " %+% paste(Q.sVars$predvars, collapse = ",") %+% ")")
    message("================================================================")
  }
  Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, 
                              predvars = Q.sVars$predvars, 
                              subset_vars = !determ.Q, 
                              estimator = getopt("Qestimator"))
  model.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData.ObsP0, savespace = TRUE)
  
  #----------------------------------------------------------------------------------
  # Create an list with model estimates, data & other information that is passed on to treatment estimation procedure
  #----------------------------------------------------------------------------------
  estinfo_list <- list(
    data = data, 
    nodes = nodes,
    TMLE.targetStep = TMLE.targetStep,
    obs.wts = obs.wts, 
    lbound = lbound,
    merged.form = merged.form, 
    model.Q.init = model.Q.init,
    Q.sVars = Q.sVars,
    f.g0 = f.g0,
    h.g0.sVars = h.g0.sVars,
    h.gstar.sVars = h.gstar.sVars,
    h.g0_GenericModel = h.g0_GenericModel,
    h.gstar_GenericModel = h.gstar_GenericModel,
    savetime.fit.hbars = savetime.fit.hbars,
    n_MCsims = n_MCsims
  ) 
  estinfo_list_g1 <- append(estinfo_list, list(f.gstar = f_gstar1))
  if (!is.null(f_gstar2)) { estinfo_list_g2 <- append(estinfo_list, list(f.gstar = f_gstar2))}
  
  #----------------------------------------------------------------------------------
  # Running MC evaluation for substitution TMLE estsimators
  #----------------------------------------------------------------------------------
  # Incl. estimate treatment mechanism f(a|E, W)) and clever covariates & targeting step
  tmle_gstar1_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g1)
  if (!is.null(f_gstar2)) {
    tmle_gstar2_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g2)
  } else {
    tmle_gstar2_out <- NULL
  }
  
  #----------------------------------------------------------------------------------
  # Create output list (estimates, as. variances, CIs)
  #----------------------------------------------------------------------------------
  EY_gstar1 <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, tmle_g_out = tmle_gstar1_out)
  EY_gstar2 <- NULL
  ATE <- NULL	
  otherInfo2 <- NULL
  if (!is.null(f_gstar2)) {
    EY_gstar2 <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, tmle_g_out = tmle_gstar2_out)
    ATE <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, 
                          tmle_g_out = tmle_gstar1_out, tmle_g2_out = tmle_gstar2_out)
  }
  message("######################################################################################")
  message("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.")
  message("######################################################################################")
  
  tmleCommunity.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)#, otherInfo1 = otherInfo1, otherInfo2 = otherInfo2)
  class(tmleCommunity.res) <- c(class(tmleCommunity.res), "tmleCommunity")
  return(tmleCommunity.res)
}
© 2017 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
