#***************************************************************************************
# Example 1: Estimate a outcome regression directly through BinaryOutModel
data(indSample.iid.bA.bY.rareJ2_list)
indSample.iid.bA.bY.rareJ2 <- indSample.iid.bA.bY.rareJ2_list$indSample.iid.bA.bY.rareJ2
N <- nrow(indSample.iid.bA.bY.rareJ2)
# speed.glm to fit regressions (it's GLMs to medium-large datasets)
tmleCom_Options(Qestimator = "speedglm__glm", maxNperBin = N)
options(tmleCommunity.verbose = TRUE)  # Print status messages 
#***************************************************************************************

#***************************************************************************************
# 1.1 Specifying outcome and predictor variables for outcome mechanism
#***************************************************************************************
# Y depends on all its parent nodes (A, W1, W2, W3, W4) 
Qform.all <- Y ~ W1 + W2 + W3 + W4 + A
Q.sVars1 <- tmleCommunity:::define_regform(regform = Qform.all)

# Equivalent way to define Q.sVars: use Anodes.lst (outcomes) & Wnodes.lst (predictors)
# node can only contain one or more of Ynode, Anodes, WEnodes, communityID and Crossnodes
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Q.sVars2 <- tmleCommunity:::define_regform(regform = NULL, Anodes.lst = nodes$Ynode, 
                                           Wnodes.lst = nodes[c("Anodes", "WEnodes")])

# Also allows to include interaction terms in regression formula  (Correct Qform)
Qform.interact <- Y ~ W1 + W2*A + W3 + W4
Q.sVars3 <- tmleCommunity:::define_regform(regform = Qform.interact)

# Alternative way to define Qform.interact 
Qform.interact2 <- Y ~ W1 + W2 + W3 + W4 + A + W2:A
Q.sVars4 <- tmleCommunity:::define_regform(regform = Qform.interact2)

#***************************************************************************************
# 1.2 Fit and predict a regression model for outcome mechanism Qbar(A, W)
#***************************************************************************************
# Create a new object of DatKeepClass that can store and munipulate the input data
OData_R6 <- DatKeepClass$new(Odata = indSample.iid.bA.bY.rareJ2, 
                             nodes = nodes, norm.c.sVars = FALSE)
# Add a vector of observation (sampling) weights that encodes knowledge of rare outcome
OData_R6$addObsWeights(obs.wts = indSample.iid.bA.bY.rareJ2_list$obs.wt.J2)

# Create a new object of RegressionClass that defines regression models
# using misspecified Qform (without interaction term) 
Qreg <- RegressionClass$new(outvar = Q.sVars1$outvars, predvars = Q.sVars1$predvars, 
                            subset_vars = (!rep_len(FALSE, N)))

# Set savespace=FALSE to save all productions during fitting, including models and data
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(data = OData_R6, savespace = FALSE)
length(m.Q.init$getY)  # 3000, the outcomes haven't been erased since savespace = FALSE
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
mean(m.Q.init$getprobA1)  # 0.02175083, bad estimate since misspecified Qform

#***************************************************************************************
# 1.3 Same as above but using Super Learner (data-adaptive algorithms)
#***************************************************************************************
# Specifying the SuperLearner library in tmleCom_Options() 
library(SuperLearner)
tmleCom_Options(SL.library = c("SL.glm", "SL.randomForest"), maxNperBin = N)
# Instead of reinitiating a RegressionClass object, change estimator directly in Qreg 
# so don't need to redefine Qestimator in tmleCom_Options()
Qreg$estimator <- "SuperLearner"

set.seed(12345)
m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(data = OData_R6, savespace = TRUE)
m.Q.init$predict(newdata = OData_R6, savespace = TRUE)
mean(m.Q.init$getprobA1)
