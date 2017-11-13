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

test_that("Assign/ return the class type of a variable.", {
  # Return the class types of variables
  expect_equal(OData_R6$get.sVar.type("A"), "contin")
  expect_equal(OData_R6$get.sVar.type("W1"), "binary")
  expect_true(all(unique(unlist(OData_R6$get.sVar.type())) %in% c("binary", "categor", "contin")))
  
  # Assign a new class type to one variable that belongs to the input data
  OData_R6.copy <- OData_R6
  OData_R6.copy$set.sVar.type(name.sVar = "W1", new.type = "contin")
  expect_equal(OData_R6.copy$get.sVar.type("W1"), "contin")
})
