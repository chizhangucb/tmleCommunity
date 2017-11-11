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

test_that("1. Invalid term name in regression formula", {
  Qform.bad <- "blabla ~ W1 + A"
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform.bad, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, 
                                obs.wts, community.wts, f_gstar1, f_gstar2),
    "Invalid term name in regression formula for 'Qform' \"blabla ~ W1 ")
  # Coudn't figure out why \"blabla ~ W1 + A\" doesn't work, so just leave it here.
  
  hform.gstar.bad <- "A ~ E1 + W5"
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform.bad, 
                                hform.g0, hform.gstar.bad, fluctuation, Qbounds, 
                                obs.wts, community.wts, f_gstar1, f_gstar2),
    "Invalid term name in regression formula for 'hform.gstar'")
})

test_that("2. The input data must be a data frame", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = as.matrix(indSample.iid.bA.bY.rareJ1), nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, 
                                community.wts, f_gstar1, f_gstar2),
    "The input data must be a data frame")
})

test_that("3. No factor column(s) allowed in the input data", {
  indSample.iid.bA.bY.rareJ1$W3 <- as.factor(indSample.iid.bA.bY.rareJ1$W3)
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, 
                                community.wts, f_gstar1, f_gstar2),
    "No factor column\\(s\\) allowed in the input data, " %+% 
      "consider removing or recoding such column\\(s\\) as strings: W3")
})

test_that("4. obs.wts must contain the same number of non-negative obs as data does", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, 
                                obs.wts = 1, community.wts, f_gstar1, f_gstar2),
    "'obs.wts', must contain the same number of non-negative observations as 'data' does")
})

test_that("5. f_gstar must contain a length (number of rows) 1 or NROW(data)", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts,
                                community.wts, f_gstar1 = 1:10, f_gstar2),
    "If 'f_gstar1' is a vector/matrix/data.frame, it should have a " %+% 
      "length \\(number of rows\\) 1 or NROW\\(data\\)")
})

test_that("6. Qbounds should have length 2", {
  expect_warning(
    tmleCommunity:::CheckInputs(data = indSample.iid.bA.bY.rareJ1, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds = 0.05, 
                                obs.wts, community.wts, f_gstar1, f_gstar2),
    "Qbounds should have length 2")
})

data("comSample.wmT.bA.bY_list", package = "tmleCommunity")
comSample.wmT.bA.bY <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"))
Qform <- "Y ~ E1 + E2 + W2 + W3 + A"
hform.g0 <- hform.gstar <- "A ~ E1 + E2 + W1" 
obs.wts <- rep(1, NROW(comSample.wmT.bA.bY))
community.wts <- 
  as.data.frame(matrix(0L, nrow = length(unique(comSample.wmT.bA.bY[, "id"])), ncol = 2))
community.wts[, 2] <-  as.vector(table(comSample.wmT.bA.bY[, "id"]))

test_that("7. community.wts must contain the same number of communities & same IDs as data", {
  warns <- capture_warnings(
    tmleCommunity:::CheckInputs(data = comSample.wmT.bA.bY, nodes, Qform, 
                                hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, 
                                community.wts[1:10, ], f_gstar1, f_gstar2))
  expect_match(warns, "'community.wts', must contain the same number of non-negative" %+% 
                 " communities as 'data' does, and has 2 columns", all = FALSE)
  expect_match(warns, "'community.wts', must contain the same community IDs as data" %+% 
                 " does \\(diff order is ok but no duplicate allows\\)", all = FALSE)
})
