context("Test DatKeepClass")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
tmleCom_Options(maxNperBin = nrow(indSample.iid.cA.cY))
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))

test_that("nodes' name can only be one or more of Ynode, Anodes, WEnodes, communityID and Crossnodes", {
  nodes.bad <- list(Ynode = "Y", Anodes = "A", badnodes = c("W1", "W2", "W3", "W4"))
  expect_error(
    expect_message(DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes.bad), "Don't recognize badnodes"), 
    "It should be a list & its names can only be one or more of Ynode, Anodes, WEnodes, communityID and Crossnodes.")
})

test_that("The length of observation weights should be the same as nrow(data)", {
  expect_error(OData_R6$addObsWeights(1:10),
    "The length of observation weights should be the same as nrow\\(data\\)")
})
