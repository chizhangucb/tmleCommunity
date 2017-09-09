### Load data
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontABinY)
dat_iidcontABinY <- sampleDat_iidcontABinY$dat_iidcontABinY
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, Crossnodes = NULL)
gvars$verbose <- TRUE  # Print status messages

# First, try to use Speedglm to fit regression (GLMs to medium-large datasets)
tmleCom_Options(Qestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
OData <- DatKeepClass$new(Odata = dat_iidcontABinY, nodes = nodes, norm.c.sVars = FALSE)
Q.sVars <- define_regform(NULL, Anodes.lst = nodes$Ynode, 
                          Wnodes.lst = nodes[c("Anodes", "Wnodes", "Enodes")])
Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars, 
                            subset_vars = (!rep_len(FALSE, nrow(dat_iidcontABinY))))
# By setting savespace = FALSE, save all data in fit step
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData, savespace = FALSE)
m.Q.init$is.fitted  # TRUE
length(m.Q.init$getY)  # 10000, the outcome observations haven't been erased since savespace = FALSE
m.Q.init$getfit$coef  # Provide cofficients from the fitting regression
# Here, wipe out any traces of saved data in predict step by setting savespace = TRUE
m.Q.init$predict(newdata = OData, savespace = TRUE)
is.null(m.Q.init$getXmat)  # TRUE, the matrix of covariates has been erased to save space
mean(m.Q.init$getprobA1)  # 0.2976

# Now, use Super Learner (data-adaptive algorithms) to fit and predict regression models
library(SuperLearner)
Qreg$estimator <- "SuperLearner"
gvars$opts$g.SL.library <- c("SL.glm", "SL.step", "SL.glm.interaction", "SL.randomForest")
set.seed(1)
# To save (RAM) space, wipe out all internal data when doing many stacked fitting regressions
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData, savespace = TRUE)
m.Q.init$predict(newdata = OData, savespace = TRUE)
mean(m.Q.init$getprobA1)  # 0.2976
