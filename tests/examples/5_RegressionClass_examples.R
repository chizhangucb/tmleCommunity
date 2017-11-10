data(indSample.iid.cA.cY_list)
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
tmleCom_Options(maxNperBin = nrow(indSample.iid.cA.cY))

OData.g0 <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)
h.g0.sVars <- tmleCommunity:::define_regform(A ~ W1 + W2 + W3 + W4) 
A.nms.g0 <- h.g0.sVars$outvars
regclass.g0 <- RegressionClass$new(outvar = A.nms.g0, predvars = h.g0.sVars$predvars,
                                   subset_vars = lapply(A.nms.g0, function(var) {var}),
                                   outvar.class = OData.g0$type.sVar[A.nms.g0])
regclass.g0$estimator  # "speedglm__glm"
regclass.g0$pool_cont  # FALSE (Don't pool across all bins of a continuous exposure)
regclass.g0$parfit     # FALSE (Don't preform parallel computing in estimation)

# Clone the parent regclass.g0 and reset to a single univariate k_i regression 
# for outcome regclass.g0$outvar[[k_i]]
k_i <- 1
reg_i <- regclass.g0$clone()
reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
genericmodels.g0.A1 <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData.g0)
