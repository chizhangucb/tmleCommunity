context("Test GenericModel")

`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
N <- nrow(indSample.iid.cA.cY)
tmleCom_Options(maxNperBin = N, nbins = 5, bin.method = "equal.mass")
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
h.g0.sVars <- tmleCommunity:::define_regform(A ~ W1 + W2 + W3 + W4)
subsets_expr <- lapply(h.g0.sVars$outvars, function(var) {var})
OData.g0 <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)
regclass.g0 <- RegressionClass$new(outvar = h.g0.sVars$outvars,
                                   predvars = h.g0.sVars$predvars,
                                   subset_vars = subsets_expr,
                                   outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    set.seed(rndseed)
    A.mu <- 0.86 * data[,"W1"] + 0.41 * data[,"W2"] - 0.34 * data[,"W3"] + 0.93 * data[,"W4"]
    untrunc.A <- rnorm(n = nrow(data), mean = A.mu + shift.const, sd = 1)
    r.new.A <- exp(0.8 * shift.const * (untrunc.A - A.mu - shift.const / 3))
    trunc.A <- ifelse(r.new.A > trunc.const, untrunc.A - shift.const, untrunc.A)
    return(trunc.A)
  }
  return(f.gstar)
}
f.gstar.corr <- define_f.gstar(shift = 2, truncBD = 10)

test_that("Different entrances to fit regression bin(s) and get the same results", {
  # Fit GenericModel directly through the parent RegressionClass
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  genericmodels.g0.A1.B2 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[2]]
  # genericmodels.g0.A1.B2$getfit$coef
  
  # Fit ContinModel through the clone of a parent RegressionClass 
  k_i <- 1  # The first exposure node
  reg_i <- regclass.g0$clone()
  reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
  genericmodels.g0.A1.cont <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData.g0)
  genericmodels.g0.A1.cont$fit(data = OData.g0, savespace = FALSE)
  genericmodels.g0.A1.cont.B2 <- genericmodels.g0.A1.alt$getPsAsW.models()[[2]]
  # genericmodels.g0.A1.alt.B2$getfit$coef
  
  # Fit GenericModel directly from genericmodels.g0.A1.alt (The first exposure node)
  subset_eval.A1 <- tmleCommunity:::def_regs_subset(genericmodels.g0.A1.cont)
  genericmodels.g0.A1.cont.Gen <- GenericModel$new(reg = subset_eval.A1)
  genericmodels.g0.A1.cont.Gen.B2 <- genericmodels.g0.A1.cont.Gen$getPsAsW.models()[[2]]
  genericmodels.g0.A1.cont.Gen.B2$fit(data = OData.g0)
  genericmodels.g0.A1.cont.Gen.B2$getfit$coef 
  
  expect_equal(genericmodels.g0.A1.B2$getfit$coef, genericmodels.g0.A1.cont.B2$getfit$coef)
  expect_equal(genericmodels.g0.A1.B2$getfit$coef, genericmodels.g0.A1.cont.Gen.B2$getfit$coef)
})

test_that("The first bin contains all zeros and so produce tiny coefficients; " %+%  
            "The last bin should contain nothing and so produce NA coefficients", {
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  lastB.Ind <- eval(regclass.g0$nbins + 2)          
  genericmodels.g0.A1.B1 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]
  genericmodels.g0.A1.Bend <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[lastB.Ind]]          
  expect_true(all(genericmodels.g0.A1.B1$getfit$coef[-1] < 10e-8))
  expect_true(all(is.na(genericmodels.g0.A1.Bend$getfit$coef)))          
})

test_that("Binarized bins contain the (approximately) same number of obs if 'equal.mass'", {
  genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
  genericmodels.g0$fit(data = OData.g0)
  genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()[[1]]
  OData.new <- OData.g0
  OData.new$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls, 
                          nbins = genericmodels.g0.A1$nbins, bin.nms = genericmodels.g0.A1$bin_nms)
  ord.sVar <- tmleCommunity:::make.ordinal(x = OData.new$get.sVar(genericmodels.g0.A1$outvar), 
                                           intervals = genericmodels.g0.A1$intrvls)
  expect_true(all(table(ord.sVar) == N / regclass.g0$nbins))
  expect_equal(sum(ord.sVar == 1), 0)
  expect_equal(sum(ord.sVar == eval(regclass.g0$nbins + 2)), 0)
})
