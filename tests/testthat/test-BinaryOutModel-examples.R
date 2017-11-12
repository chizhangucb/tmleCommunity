context("Test BinaryOutModel using different estimators")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.bA.bY.rareJ1_list", package = "tmleCommunity")
indSample.iid.bA.bY.rareJ1 <- indSample.iid.bA.bY.rareJ1_list$indSample.iid.bA.bY.rareJ1
N <- nrow(indSample.iid.bA.bY.rareJ1)
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Q.sVars <- tmleCommunity:::define_regform(regform = Y ~ W1 + W2 + W3 + W4 + A)

test_that("Using glm__glm when setting Qestimator = 'glm__glm'", {
  tmleCom_Options(Qestimator = "glm__glm", maxNperBin = N)
  Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars,
                              subset_vars = (!rep_len(FALSE, N)))
  m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(data = OData_R6)
  m.Q.init$predict(newdata = OData_R6)
  expect_equal(m.Q.init$getfit$fitfunname, "glm")
})
