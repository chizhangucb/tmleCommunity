\dontrun{
#***************************************************************************************
# Example 1: Defining and modeling P(A | W) with continuous A
data(indSample.iid.cA.cY_list)
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
tmleCom_Options(gestimator = "speedglm__glm", maxNperBin = nrow(indSample.iid.cA.cY),
                bin.method = "equal.mass", nbins = 10)
options(tmleCommunity.verbose = FALSE)  # Don't print status messages
#***************************************************************************************

#***************************************************************************************
# 1.1 Defining new R6 objects of DatKeepClass and RegressionClass and GenericModel
#***************************************************************************************
OData.g0 <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)
h.g0.sVars <- tmleCommunity:::define_regform(NULL, Anodes.lst = nodes$Anodes, 
                                             Wnodes.lst = nodes$WEnodes)
subsets_expr <- lapply(h.g0.sVars$outvars, function(var) {var}) 
regclass.g0 <- RegressionClass$new(outvar = h.g0.sVars$outvars,
                                   predvars = h.g0.sVars$predvars,
                                   subset_vars = subsets_expr,
                                   outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)

#***************************************************************************************
# 1.2 Details regarding the GenericModel of the first exposure variable  
#***************************************************************************************
genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1` 
genericmodels.g0$getPsAsW.models()$`P(A|W).2`  # NULL as only one A in the input data
# Defining the number and positions of the bins via arguments in tmleCom_Options()
genericmodels.g0.A1$intrvls  
# Creating a matrix of dummy bin indicators for continuous A
OData.g0$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, 
                       intervals = genericmodels.g0.A1$intrvls, 
                       nbins = genericmodels.g0.A1$reg$nbins, 
                       bin.nms = genericmodels.g0.A1$bin_nms)  
bin.ind.mat <- OData.g0$dat.bin.sVar
colSums(bin.ind.mat, na.rm = TRUE)   # Each bin has 1000 obs as "equal.mass" with 10 bins

#***************************************************************************************
# 1.3 Fitting regression models for the first exposure variable  
#***************************************************************************************
genericmodels.g0.A1$fit(data = OData.g0)
genericmodels.g0.A1.B2 <- genericmodels.g0.A1$getPsAsW.models()$`P(A|W).2`  # 2nd bin
genericmodels.g0.A1.B2$getfit$coef
}
