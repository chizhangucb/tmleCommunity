context("Test BinaryOutModel using different estimators")
context("Test for saving RAM during fitting and predicting process")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
N <- nrow(indSample.iid.cA.cY)
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Q.sVars <- tmleCommunity:::define_regform(regform = Y ~ W1 + W2 + W3 + W4 + A)
h.g0.sVars <- tmleCommunity:::define_regform(A ~ W1 + W2 + W3 + W4)
subsets_expr <- lapply(h.g0.sVars$outvars, function(var) {var})
OData.g0 <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)

#**********************************************  
# Test 1 Different estimation algorithms
#********************************************** 
test_that("Using glm when setting Qestimator = 'glm__glm'", {
  # Use glm without pooling of bins
  tmleCom_Options(Qestimator = "glm__glm", maxNperBin = N)
  regclass.g0 <- RegressionClass$new(
    outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
    subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1`
  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
    expect_equal(genericmodels.g0.A1$getPsAsW.models()[[i]]$estimator, "glm__glm")
  }
  
  # Use glm with pooling of bins (glm.long)
  # tmleCom_Options(Qestimator = "glm__glm", maxNperBin = N, poolContinVar = TRUE)
  # regclass.g0 <- RegressionClass$new(
  #   outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
  #   subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
  # genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  # genericmodels.g0$fit(data = OData.g0)
  # h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
})

test_that("Using speedglm when setting Qestimator = 'speedglm__glm'", {
  tmleCom_Options(Qestimator = "speedglm__glm", maxNperBin = N)
  regclass.g0 <- RegressionClass$new(
    outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
    subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1`
  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
    expect_equal(genericmodels.g0.A1$getPsAsW.models()[[i]]$estimator, "speedglm__glm")
  }
})

test_that("Using SuperLearner when setting Qestimator = 'SuperLearner'", {
  require("SuperLearner")
  tmleCom_Options(Qestimator = "SuperLearner", maxNperBin = N, 
                  SL.library = c("SL.glm", "SL.stepAIC", "SL.bayesglm"))
  regclass.g0 <- RegressionClass$new(
    outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
    subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1`
  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
    expect_equal(genericmodels.g0.A1$getPsAsW.models()[[i]]$estimator, "SuperLearner")
  }
  # For the first and last bin, no much obs in it, so SL fails & downgrade to speedglm 
  expect_equal(genericmodels.g0.A1$getPsAsW.models()[[1]]$getfit$fitfunname, "speedglm")
})

#test_that("Using h2o & h2oEnsemble when setting gestimator = 'h2o__ensemble'", {
#  require("h2o"); require("h2oEnsemble")
#  tmleCom_Options(Qestimator = "h2o__ensemble", maxNperBin = N, 
#                  h2ometalearner = "h2o.glm.wrapper",
#                  h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"))
#  regclass.g0 <- RegressionClass$new(
#    outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
#    subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
#  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
#  genericmodels.g0$fit(data = OData.g0)
#  h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
#  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1`
#  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
#    expect_equal(genericmodels.g0.A1$getPsAsW.models()[[i]]$estimator, "h2o__ensemble")
#  }
#  # For the last (& first) bin, no much obs in it, so h2o fails & downgrade to speedglm 
#  expect_equal(genericmodels.g0.A1$getPsAsW.models()[[7]]$getfit$fitfunname, "speedglm")
#})

#**********************************************  
# Test 2 Test for saving RAM 
#********************************************** 
OData_R6 <- DatKeepClass$new(Odata = subset(indSample.iid.cA.cY, select=-Y),
                             nodes = nodes[c("Anodes", "WEnodes")], norm.c.sVars = FALSE)
OData_R6$nodes <- nodes
obsYvals <- indSample.iid.cA.cY[, nodes$Ynode]
ab <- range(obsYvals, na.rm=TRUE)
obsYvals.bd <- (obsYvals-ab[1]) / diff(ab)
OData_R6$addYnode(YnodeVals = obsYvals.bd, det.Y = FALSE)
Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars, 
                            subset_vars = (!rep_len(FALSE, OData_R6$nobs)))

# Test 2.1 Wipe out any traces of saved data in both fit and predict steps
test_that("Wiping out any traces of saved data after fitting & predicting regression", {
  m.Q.init_wipefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData_R6, savespace = TRUE)
  m.Q.init_wipefit_wipePred <- m.Q.init_wipefit$clone(deep = T)
  m.Q.init_wipefit_wipePred$predict(newdata = OData_R6, savespace = TRUE)
  expect_true(is.null(c(m.Q.init_wipefit$getXmat, m.Q.init_wipefit$getY)))
  expect_true(is.null(m.Q.init_wipefit_wipePred$getXmat))
})

# Test 2.2 Save all data in fit step but wipe out any traces of saved data in predict step
test_that("saving all after fitting regression but wiping out any traces of saved data after prediction", {
  m.Q.init_savefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData_R6, savespace = FALSE)
  m.Q.init_savefit_wipePred <- m.Q.init_savefit$clone(deep = T)
  m.Q.init_savefit_wipePred$predict(newdata = OData_R6, savespace = TRUE)
  expect_length(m.Q.init_savefit$getY, 10000)
  expect_true(class(m.Q.init_savefit$getXmat) == "matrix")
  expect_true(is.null(m.Q.init_savefit_wipePred$getXmat))
})

# Test 2.3 Wipe out any traces of saved data in fit step but  Save all data in predict step
test_that("saving all generated data after fitting & predicting regression", {
  m.Q.init_savefit <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData_R6, savespace = FALSE)
  m.Q.init_savefit_savePred <- m.Q.init_savefit$clone(deep = T)
  m.Q.init_savefit_savePred$predict(newdata = OData_R6, savespace = FALSE)
  expect_equivalent(m.Q.init_savefit$getY, m.Q.init_savefit_savePred$getY)
  expect_true(class(m.Q.init_savefit_savePred$getXmat) == "matrix")
})
