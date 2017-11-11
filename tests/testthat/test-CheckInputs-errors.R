#***************************************************************************************
context("Test CheckInput Error Handling")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.bA.bY.rareJ1_list", package = "tmleCommunity")
indSample.iid.bA.bY.rareJ1 <- indSample.iid.bA.bY.rareJ1_list$indSample.iid.bA.bY.rareJ1
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
Qform <- "Y ~ W1 + W2*A + W3 + W4"
hform.g0 <- hform.gstar <- "A ~ W1 + W2 + W3 + W4"
fluctuation <- "logistic"
Qbounds <- NULL
obs.wts <- rep(1, NROW(indSample.iid.bA.bY.rareJ1))
community.wts <- NULL
f_gstar1 <- NULL
f_gstar2 <- NULL

test_that("Invalid term name in regression formula", {
  Qform.bad <- "blabla ~ W1 + A"
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform.bad, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, 
                                obs.wts, community.wts, f_gstar1, f_gstar2),
    "Invalid term name in regression formula for 'Qform'")
  
  hform.gstar.bad <- "A ~ E1 + W5"
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform.bad, 
                                hform.g0, hform.gstar.bad, fluctuation, Qbounds, 
                                obs.wts, community.wts, f_gstar1, f_gstar2),
    "Invalid term name in regression formula for 'hform.gstar'")
})

test_that("The input data must be a data frame", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = as.matrix(indSample.iid.bA.bY.rareJ1), nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, 
                                community.wts, f_gstar1, f_gstar2),
    "The input data must be a data frame")
})

test_that("No factor column(s) allowed in the input data", {
  indSample.iid.bA.bY.rareJ1$W3 <- as.factor(indSample.iid.bA.bY.rareJ1$W3)
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, 
                                community.wts, f_gstar1, f_gstar2),
    "No factor column(s) allowed in the input data")
})

test_that("obs.wts must contain the same number of non-negative obs as data does", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, 
                                obs.wts = 1, community.wts, f_gstar1, f_gstar2),
    "'obs.wts', must contain the same number of non-negative observations as 'data' does")
})

test_that("f_gstar must contain a length (number of rows) 1 or NROW(data)", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts,
                                community.wts, f_gstar1 = 1:10, f_gstar2),
    "If 'f_gstar1' is a vector/matrix/data.frame, it should have a " %+% 
      "length \\(number of rows\\) 1 or NROW\\(data\\)")
})
