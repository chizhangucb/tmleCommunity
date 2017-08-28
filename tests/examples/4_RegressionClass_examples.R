data(sampleDat_iidcontAContY)
dat_iidcontAContY <- sampleDat_iidcontAContY$dat_iidcontAContY
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, Crossnodes = NULL) 
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", parfit = FALSE,
                maxNperBin = nrow(dat_iidcontAContY), poolContinVar = FALSE)
OData <- DatKeepClass$new(Odata = dat_iidcontAContY, nodes = nodes, norm.c.sVars = FALSE)

# Parse the formulae for actual covariate names in A, W & E
# If the formula is not specified (i.e., NULL), variable(s) input into Anodes.lst will be treated as outcome
# variable(s), and all variables input into Wnodes.lst will be used as explanatory variable(s).
g.sVars.null <- define_regform(NULL, Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[c("Wnodes", "Enodes")])
g.sVars.null$predvars  # "W1" "W2" "W3" "W4"

# When use correctly specified g mechanism
gform.corr <- "A ~ W1 + W2 + W3 + W4" 
g.sVars.corr <- define_regform(as.formula(gform.corr), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[3:4]) 
g.sVars.corr$predvars  # "W1" "W2" "W3" "W4"

# When use incorrectly specified g mechanism
gform.mis <- "A ~ W2 * W3 * W4"
g.sVars.mis <- define_regform(as.formula(gform.mis), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[3:4])
g.sVars.mis$predvars  # "W2" "W3" "W4" "W2:W3" "W2:W4" "W3:W4" "W2:W3:W4"

A_nms_g0 <- g.sVars.corr$outvars
regclass.g0 <- RegressionClass$new(outvar = A_nms_g0, predvars = g.sVars.corr$predvars,
                                   subset_vars = lapply(A_nms_g0, function(var) {var}) ,
                                   outvar.class = OData$type.sVar[A_nms_g0])
regclass.g0$estimator  # "speedglm__glm"
regclass.g0$pool_cont  # FALSE (Don't pool across all bins of a continuous exposure)
regclass.g0$parfit     # FALSE (Don't preform parallel computing in estimation)

# Use a RegressionClass object to construct a GenericModel object
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData)
genericmodels.g0$is.fitted  # FALSE

# Details regarding A1 (i.e., the first exposure)
genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1` # genericmodels.g0$getPsAsW.models()[[1]]
genericmodels.g0.A1$intrvls
# Create a matrix of dummy bin indicators for continuous sVar 
OData$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls, 
                    nbins = genericmodels.g0.A1$reg$nbins, bin.nms = genericmodels.g0.A1$bin_nms)  
bin.ind.mat <- OData$dat.bin.sVar
colSums(bin.ind.mat, na.rm = T)
#  A_B.1  A_B.2  A_B.3  A_B.4  A_B.5  A_B.6  A_B.7  A_B.8  A_B.9 A_B.10 A_B.11 A_B.12 
#    0   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000      0 
genericmodels.g0$fit(data = OData)
genericmodels.g0.A1.B1 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]
genericmodels.g0.A1.B2 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[2]]
genericmodels.g0.A1.B2$getfit$coef  # the same as genericmodels.g0.A1.alt.B2$getfit$coef
#  Intercept         W1         W2         W3         W4 
# -0.6758011 -1.6536866 -0.7941279  0.5959984 -1.6902898

# Set regclass.g0 to a single univariate k_i regression for outcome regclass.g0$outvar[[k_i]]
k_i <- 1
reg_i <- regclass.g0$clone()
reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
genericmodels.g0.A1.alt <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData)
