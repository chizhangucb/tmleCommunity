fitGenericDensity <- function(data, Anodes, Wnodes, gform, f_gstar, h.gstar_GenericModel = NULL, 
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
  nobs <- OData.ObsP0$nobs
  
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
