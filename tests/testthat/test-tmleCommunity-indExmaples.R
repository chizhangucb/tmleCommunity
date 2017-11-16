context("Test tmleCommunity on individual-level exposure")

# ---------------------------------------------------------------------------------
# TEST SET 1. TESTS FOR FITTING BINARY EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
get.iid.dat.Abin <- function(ndata = 100000, rndseed = NULL, is.Y.bin = TRUE) {
  require(simcausal)
  D <- DAG.empty()
  D <- D + 
    node("W1", distr = "rbern", prob = 0.5) +
    node("W2", distr = "rbern", prob = 0.3) +
    node("W3", distr = "rnorm", mean = 0, sd = 0.3) +
    node("W4", distr = "runif", min = 0, max = 1) +
    node("W3W4", distr = "rconst", const = W3 * W4) +
    node("A", distr = "rbern", prob = plogis(0.86 * W1 + W2 + 1.32 * W3 - W4 - 0.45 *W3W4 + 0.1)) +
    node("Y", distr = "rnorm", mean = (2.8 * A + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) + 
    node("Y1", distr = "rnorm", mean = (2.8 * 1 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) +
    node("Y0", distr = "rnorm", mean = (2.8 * 0 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1)
  
  D <- set.DAG(D)
  Odata <- sim(D, n = ndata, rndseed = rndseed)
  psi0.Y <- mean(Odata$Y1) - mean(Odata$Y0)
  print("psi0.Y: " %+% psi0.Y)  
  return(list(psi0.Y = psi0.Y, Odata = Odata))
}

`%+%` <- function(a, b) paste0(a, b)
ndata <- 10000
rndseed <- 12345
indSample.iid.bA.cY <- get.iid.dat.Abin(ndata = ndata, rndseed = rndseed)$Odata
N <- NROW(indSample.iid.bA.cY)
indSample.iid.bA.cY <- indSample.iid.bA.cY[, c("W1", "W2", "W3", "W4", "A", "Y")]

Qform.corr <- "Y ~ W1 + W2 + W3 + W4 + A" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 * W4"  # correct g
gform.mis <- "A ~ W1 + W3"  # incorrect g

#******************************
# Test 1.1 speedglm & glm 
#******************************
## Test 1.1.1 Correctly specified Qform & gform (+ tmle.intercept + equal.mass + 5 nbins)
test_that("fit TMLE for binary A with speedglm, with correctly specified Qform & gform", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  variances <- tmleCom_res$ATE$vars
  expect_equal(estimates["tmle", ], 2.799876, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 2.783021, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 2.779068, tolerance = 1e-6) 
  expect_equal(as.vector(variances), c(0.000491327, 0.027397109, 0.000491327))
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "environment") 
})

## Test 1.1.2 Misspecified Qform (+ tmle.intercept + equal.mass + 5 nbins)
test_that("fit TMLE for binary A with speedglm, with misspecified Qform", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.mis, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  expect_equal(estimates["tmle", ], 2.801727, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 2.783021, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 3.696836, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["gcomp", ] - psi0.Y))
})

## Test 1.1.3 Misspecified gform (+ tmle.intercept + equal.mass + 5 nbins)
test_that("fit TMLE for binary A with speedglm, with misspecified gform", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  expect_equal(estimates["tmle", ], 2.777445, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 3.307761, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 2.779068, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["iptw", ] - psi0.Y))
})

## Test 1.1.4 Misspecified gform (+ tmle.intercept + equal.mass + 20 nbins)
test_that("fit TMLE for binary A with speedglm, with misspecified gform (20 bins)", {
  tmleCom_Options(maxNperBin = N, nbins = 20)
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  expect_equal(estimates["tmle", ], 2.777445, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 3.307761, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 2.779068, tolerance = 1e-6) 
})

#******************************
## Test 1.2 SuperLearner
#******************************
## Test 1.2.1 SuperLearners with only main terms (+ tmle.intercept + equal.mass + 5 nbins)
test_that("fit TMLE for binary A with SuperLearner, using SL.glm, SL.step, SL.glm.interaction,
          when both Qform and gform only use main terms", {
  require("SuperLearner")
  tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = N,
                  SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"))
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0,
                  Qform = NULL, hform.g0 = NULL, hform.gstar = NULL, rndseed = 12345)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  expect_equal(estimates["tmle", ], 2.799282, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 2.801273, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 2.780339, tolerance = 1e-6) 
})

#******************************
# Test 1.3 h2o
#******************************
## Test 1.3.1 h2o.glm.wrapper with only main terms (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for binary A with h2o, using h2o.glm.wrapper and main term formulae", {
  require("h2oEnsemble")           
  tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = N,
                  h2olearner = c("h2o.glm.wrapper"), nbins = 5)
  tmleCom_res <- 
    tmleCommunity(data = indSample.iid.bA.cY, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, rndseed = 12345)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.80026 
  expect_equal(estimates["tmle", ], 2.800851, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 2.874219, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 2.800498, tolerance = 1e-6) 
  h2o.removeAll()
})

# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN IID DATA
# --------------------------------------------------------------------------------- 
