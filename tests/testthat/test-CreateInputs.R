context("Test CreateInputs")

data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY

test_that("Outcome outside of [0, 1] is bounded back when maptoYstar=TRUE", {
  alpha <- 0.995
  obsYvals <- indSample.iid.cA.cY$Y
  boundedY.list <- tmleCommunity:::CreateInputs(obsYvals, Qbounds = NULL, alpha = alpha, maptoYstar = TRUE)
  expect_equal(boundedY.list$Ystar * diff(boundedY.list$ab) + boundedY.list$ab[1], obsYvals)
  expect_equal(boundedY.list$Qbounds, c(1 - alpha, alpha))
})

test_that("Outcome outside of [0, 1] is unchanged when maptoYstar=FALSE", {
  alpha <- 0.995
  obsYvals <- indSample.iid.cA.cY$Y
  boundedY.list <- tmleCommunity:::CreateInputs(obsYvals, Qbounds = NULL, alpha = alpha, maptoYstar = FALSE)
  expect_equal(boundedY.list$Ystar, obsYvals)
  expect_equal(boundedY.list$ab, c(0, 1))
})
