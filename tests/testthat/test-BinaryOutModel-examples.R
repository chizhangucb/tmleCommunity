context("Test BinaryOutModel using different estimators")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.bA.bY.rareJ1_list", package = "tmleCommunity")
indSample.iid.bA.bY.rareJ1 <- indSample.iid.bA.bY.rareJ1_list$indSample.iid.bA.bY.rareJ1
N <- nrow(indSample.iid.bA.bY.rareJ1)
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Q.sVars <- tmleCommunity:::define_regform(regform = Y ~ W1 + W2 + W3 + W4 + A)

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

test_that("Using h2o & h2oEnsemble when setting gestimator = 'h2o__ensemble'", {
  require("h2o"); require("h2oEnsemble")
  tmleCom_Options(Qestimator = "h2o__ensemble", maxNperBin = N, 
                  h2ometalearner = "h2o.glm.wrapper",
                  h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"))
  regclass.g0 <- RegressionClass$new(
    outvar = h.g0.sVars$outvars, predvars = h.g0.sVars$predvars,
    subset_vars = subsets_expr, outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  h_gN <- genericmodels.g0$predictAeqa(newdata = OData.g0)
  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1`
  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
    expect_equal(genericmodels.g0.A1$getPsAsW.models()[[i]]$estimator, "h2o__ensemble")
  }
  # For the first (& last) bin, no much obs in it, so h2o fails & downgrade to speedglm 
  expect_equal(genericmodels.g0.A1$getPsAsW.models()[[1]]$getfit$fitfunname, "speedglm")
})
