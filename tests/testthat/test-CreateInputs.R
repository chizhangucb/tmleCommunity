context("Test CreateInputs")

data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY

test_that("Continuous outcome will always be bounded into [0, 1] when maptoYstar=TRUE", {
  obsYvals <- indSample.iid.cA.cY$Y
  boundedY.list <- CreateInputs(obsYvals, Qbounds = NULL, alpha = 0.995, maptoYstar = TRUE)
  expect_equal(boundedY.list$Ystar * diff(boundedY.list$ab) + boundedY.list$ab[1], obsYvals)
})
