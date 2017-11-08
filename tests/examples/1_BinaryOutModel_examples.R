#***************************************************************************************
# Example 1: Estimate a outcome regression directly through BinaryOutModel
data(indSample.iid.cA.bY_list)
indSample.iid.cA.bY <- indSample.iid.cA.bY_list$indSample.iid.cA.bY
# speed.glm to fit regressions (it's GLMs to medium-large datasets)
tmleCom_Options(Qestimator = "speedglm__glm", maxNperBin = nrow(indSample.iid.cA.bY))
gvars$verbose <- TRUE  # Print status messages (global setting)
#***************************************************************************************

#***************************************************************************************
# 1.1 Specifying outcome and predictor variables for outcome mechanism
#***************************************************************************************
# Y depends on all its parent nodes (A, W1, W2, W3, W4)
Qform.corr <- Y ~ W1 + W2 + W3 + W4 + A
# node can only contain one or more of Ynode, Anodes, WEnodes, communityID and Crossnodes
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Q.sVars1 <- define_regform(regform = Qform.corr)

# Equivalent way to define Q.sVars is to provide Anodes.lst (outcomes) & Wnodes.lst 
# (predictors) (No need for regform)
Q.sVars2 <- define_regform(regform = NULL, Anodes.lst = nodes$Ynode, 
                           Wnodes.lst = nodes[c("Anodes", "WEnodes")])

# Also allows to include interaction terms in regression formula
Qform.interact <- Y ~ W1 + W3 * A
Q.sVars4 <- define_regform(regform = Qform.interact)

# Alternative way to define Qform.interact 
Qform.interact2 <- Y ~ W1 + W3 + A + W3:A
Q.sVars5 <- define_regform(regform = Qform.interact2)

#***************************************************************************************
# 1.2 Fit and predict a regression model for outcome mechanism Qbar(A, W)
#***************************************************************************************
# Create a new object of DatKeepClass that can store and munipulate the input data
OData_R6 <- DatKeepClass$new(Odata = indSample.iid.cA.bY, 
                             nodes = nodes, norm.c.sVars = FALSE)
# Create a new object of RegressionClass that defines regression models
Qreg <- RegressionClass$new(outvar = Q.sVars1$outvars, predvars = Q.sVars1$predvars, 
                            subset_vars = (!rep_len(FALSE, nrow(indSample.iid.cA.bY))))

# Set savespace=FALSE to save all productions during fitting, including models and data
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(data = OData_R6, savespace = FALSE)
length(m.Q.init$getY)  # 10000, the outcomes haven't been erased since savespace = FALSE
head(m.Q.init$getXmat)  # the predictor matrix is kept since savespace = FALSE
m.Q.init$getfit$coef  # Provide cofficients from the fitting regression
m.Q.init$is.fitted  # TRUE

# Now fit the same regression model but set savespace to TRUE (only fitted model left)
# Need to set overwrite to TRUE to avoid error when m.Q.init is already fitted
m.Q.init <- m.Q.init$fit(overwrite = TRUE, data = OData_R6, savespace = TRUE)
all(is.null(m.Q.init$getXmat), is.null(m.Q.init$getY))  # TRUE, all wiped out

# Set savespace = TRUE to wipe out any traces of saved data in predict step
m.Q.init$predict(newdata = OData_R6, savespace = TRUE)
is.null(m.Q.init$getXmat)  # TRUE, the covariates matrix has been erased to save RAM space
mean(m.Q.init$getprobA1)  # 0.3351

#***************************************************************************************
# 1.3 Same as above but using Super Learner (data-adaptive algorithms)
#***************************************************************************************
# Instead of starting over from tmleCom_Options(), there is a more direct way to 
# specify the estimator used in Qreg and SuperLearner library
Qreg$estimator <- "SuperLearner"
# Another way to specify the prediction algorithms in SL library
gvars$opts$g.SL.library <- c("SL.glm", "SL.glm.interaction", "SL.randomForest")

set.seed(12345)
library(SuperLearner)
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(data = OData_R6, savespace = TRUE)
m.Q.init$predict(newdata = OData_R6, savespace = TRUE)
mean(m.Q.init$getprobA1)  # 0.3351
