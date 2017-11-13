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

test_that("Assign/ return/ check the class type of a variable.", {
  # Return the class types of variables
  expect_equal(OData_R6$get.sVar.type("A"), "contin")
  expect_equal(OData_R6$get.sVar.type("W1"), "binary")
  expect_true(all(unique(unlist(OData_R6$get.sVar.type())) %in% c("binary", "categor", "contin")))
  
  # Assign a new class type to one variable that belongs to the input data
  OData_R6.copy <- OData_R6
  OData_R6.copy$set.sVar.type(name.sVar = "W1", new.type = "contin")
  expect_equal(OData_R6.copy$get.sVar.type("W1"), "contin")
    
  # Check the class types of variables
  expect_true(OData_R6$is.sVar.cont("A"))
  expect_true(OData_R6$is.sVar.bin("W1"))
  expect_false(OData_R6$is.sVar.cat("W3"))
})

test_that("Three binning methods for continuous/ categorical sVar", {
  OData_R6.copy <- OData_R6
  
  # "equal.mass" will create bins with (approximately) the same number of obs
  intrvls <- OData_R6.copy$detect.sVar.intrvls(
    name.sVar = "A", nbins = 5, bin_bymass = TRUE, bin_bydhist = FALSE, max_nperbin = N)
  OData_R6.copy$binirize.sVar(name.sVar = "A", intervals = intrvls, nbins = eval(5 + 2), 
    bin.nms = OData_R6.copy$bin.nms.sVar("A", eval(5 + 2)))
  expect_equal(
    as.vector(colSums(OData_R6.copy$dat.bin.sVar, na.rm = TRUE)[2:6]), rep(N/5, 5))
  
  # "equal.len" will create bins with the same length
  intrvls <- OData_R6.copy$detect.sVar.intrvls(
    name.sVar = "A", nbins = 5, bin_bymass = FALSE, bin_bydhist = FALSE, max_nperbin = N)
  expect_true(isTRUE(all.equal(max(diff(intrvls)[2:5]), min(diff(intrvls)[2:5]))))
  
  # "dhist" is a compromise between "equal.mass" and "equal.len"
  intrvls <- OData_R6.copy$detect.sVar.intrvls(
    name.sVar = "A", nbins = 5, bin_bymass = FALSE, bin_bydhist = TRUE, max_nperbin = N)
  OData_R6.copy$binirize.sVar(name.sVar = "A", intervals = intrvls, nbins = eval(5 + 2), 
                              bin.nms = OData_R6.copy$bin.nms.sVar("A", eval(5 + 2)))
  expect_equal(sum(OData_R6.copy$dat.bin.sVar, na.rm = TRUE), N)
})

test_that("Making a new dataframe under stochastic intervention f.g_fun", {
  OData_R6.copy <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)
  define_f.gstar <- function(shift.rate) {
    eval(shift.rate)
    f.gstar <- function(data, ...) { data[, "A"] * shift.rate }
    return(f.gstar)
  }
  OData_R6.copy$make.dat.sVar(p = 1, f.g_fun = define_f.gstar(0.4))
  # Baseline covariates doesn't change
  expect_equal(OData_R6.copy$get.dat.sVar(covars = nodes$WEnodes), 
               OData_R6$get.dat.sVar(covars = nodes$WEnodes))
  # Exposure values change via stochastic intervention function
  expect_equal(OData_R6.copy$get.sVar("A"), OData_R6$get.sVar("A") * 0.4)
    
  # If f.g_fun is a vector or dataframe/ matrix, its length (nrow) should be 1 or NROW(data)
  expect_error(OData_R6.copy$make.dat.sVar(p = 1, f.g_fun = 1:2), 
               "f_gstar1/f_gstar2 must be either a function or a vector of length nrow\\(data\\) or 1")
    
  # If f.g_fun is a function, it must contain a named argument 'data'
  expect_error(OData_R6.copy$make.dat.sVar(p = 1, f.g_fun = mean), 
               "functions f_gstar1 / f_gstar2 must have a named argument 'data'")
})
