#***************************************************************************************
# Example 1: Defining and modeling P(A | W) with continuous A
data(indSample.iid.cA.bY_list)
indSample.iid.cA.bY <- indSample.iid.cA.bY_list$indSample.iid.cA.bY
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
tmleCom_Options(gestimator = "speedglm__glm", maxNperBin = nrow(indSample.iid.cA.bY))
gvars$verbose <- TRUE  # Print status messages (global setting)
#***************************************************************************************

OData_R6 <- DatKeepClass$new(Odata = indSample.iid.cA.bY, nodes = nodes)

h.g0.sVars <- define_regform(NULL, Anodes.lst = nodes$Anodes, Wnodes.lst = nodes$WEnodes)
sA_nms_g0 <- h.g0.sVars$outvars
subsets_expr <- lapply(sA_nms_g0, function(var) {var}) 
regclass.g0 <- RegressionClass$new(outvar = h.g0.sVars$outvars,
                                   predvars = h.g0.sVars$predvars,
                                   subset_vars = subsets_expr,
                                   outvar.class = OData$type.sVar[sA_nms_g0])
