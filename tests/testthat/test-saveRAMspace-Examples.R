context("Test for saving RAM during fitting and predicting process")


# FOCUS ON using arguments in GenericModelClasses.R
# Show choices of keeping all data/models during estimation OR wiping them out
library(testthat)
gvars$verbose <- TRUE

# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontABinY)
dat_iidcontABinY <- sampleDat_iidcontABinY$dat_iidcontABinY
psi0.Y <- mean(dat_iidcontABinY$Y)  # 0.29154
psi0.Ygstar <- mean(dat_iidcontABinY$Y.gstar)  # 0.31627
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, Crossnodes = NULL)

## Create an R6 object that stores and manages the input data, later passed on to estimation algorithm(s)
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
OData <- DatKeepClass$new(Odata = dat_iidcontABinY, nodes = nodes, norm.c.sVars = FALSE)
nobs <- OData$nobs
obsYvals <- dat_iidcontABinY[, nodes$Ynodes]
OData$addYnode(YnodeVals = obsYvals)
Q.sVars <- define_regform(NULL, Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "Wnodes", "Enodes")])
Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars, subset_vars = (!rep_len(FALSE, nobs)))

######################################### 
## Test 1 Wipe out any traces of saved data in both fit and predict steps
######################################### 
test_that("wiping out any traces of saved data after fitting & predicting regression", {
  m.Q.init_wipefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData,  savespace = TRUE)
  m.Q.init_wipefit_wipePred <- m.Q.init_wipefit$clone(deep = T)
  m.Q.init_wipefit_wipePred$predict(newdata = OData, savespace = TRUE)
  expect_true(is.null(c(m.Q.init_wipefit$getXmat, m.Q.init_wipefit$getY)))
  expect_true(is.null(m.Q.init_wipefit_wipePred$getXmat))
})

######################################### 
## Test 2 Save all data in fit step but wipe out any traces of saved data in predict step
######################################### 
test_that("saving all after fitting regression but wiping out any traces of saved data after prediction", {
  m.Q.init_savefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData,  savespace = FALSE)
  m.Q.init_savefit_wipePred <- m.Q.init_savefit$clone(deep = T)
  m.Q.init_savefit_wipePred$predict(newdata = OData, savespace = TRUE)
  expect_length(m.Q.init_savefit$getY, 10000)
  expect_true(class(m.Q.init_savefit$getXmat) == "matrix")
  expect_true(is.null(m.Q.init_savefit_wipePred$getXmat))
})

######################################### 
## Test 3 Wipe out any traces of saved data in fit step but  Save all data in predict step
######################################### 
test_that("saving all generated data after fitting & predicting regression", {
  m.Q.init_savefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData,  savespace = FALSE)
  m.Q.init_savefit_savePred <- m.Q.init_savefit$clone(deep = T)
  m.Q.init_savefit_savePred$predict(newdata = OData, savespace = FALSE)
  expect_equivalent(m.Q.init_savefit$getY, m.Q.init_savefit_savePred$getY)
  expect_true(class(m.Q.init_savefit_savePred$getXmat) == "matrix")
})
