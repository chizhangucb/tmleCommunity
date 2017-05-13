fitGenericDensity <- function(data, Anodes, Wnodes, gform, f_gstar, f.g0 = NULL, gform_GenericModel = NULL, 
                              gestimator, obs.wts = NULL, rndseed = NULL, verbose = TRUE) {
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  if (!missing(gestimator)) {
    if (!(gestimator %in% gvars$opts.allowedVals[["gestimator"]])) 
      stop("gestimator must be one of: " %+% paste0(gvars$opts.allowedVals[["gestimator"]], collapse=", "))
  } else {
    gestimator <- getopt("gestimator")
  }
  gvars$verbose <- verbose
  message("Running fitGenericDensity with the following settings from tmleCom_Options(): "); str(gvars$opts)
  
  if (is.null(obs.wts)) obs.wts <- rep(1, nrow(data))
 
  ## Check if any unexpected inputs
  nodes <- list(Anodes = Anodes, Wnodes = Wnodes, Enodes = Enodes, Crossnodes = NULL)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, NULL, gform, NULL, "logistic", NULL, obs.wts)) stop()
  
  g.sVars <- define_regform(as.formula(gform), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes$Wnodes)
  
  if (verbose) {
    print("Input regression gform (P(A|W,E) under g0): " %+% gform)
    print("Derived regression gform (P(A|W,E) under g0): "); str(g.sVars)
  }
  
}
