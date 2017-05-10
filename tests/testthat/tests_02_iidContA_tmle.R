# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by binning, conditional on covariates
# ---------------------------------------------------------------------------------
library(testthat)
gvars$verbose <- TRUE

# ---------------------------------------------------------------------------------
# Test 1. The TMLE/ IPTW/ GCOMP estimator for binary Y 
# --------------------------------------------------------------------------------- 
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontABinY)
dat_iidcontABinY <- sampleDat_iidcontABinY$dat_iidcontABinY
head(dat_iidcontABinY)
psi0.Y <- sampleDat_iidcontABinY$psi0.Y  # 0.291398
psi0.Ygstar <- sampleDat_iidcontABinY$psi0.Ygstar  # 0.316274

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
f.gstar <- define_f.gstar(shift = sampleDat_iidcontABinY$shift.val, truncBD = sampleDat_iidcontABinY$truncBD, 
                          rndseed = sampleDat_iidcontABinY$rndseed)

Qform.corr <- "Y ~ W1 + W2 + W3 + W4 + A" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g
gform.mis <- "A ~ W2 + W3"  # correct g

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
## Test 1.1.1 Correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when Qform & gform are correctly specified", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3280404, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.3278714, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3293320, tolerance = 1e-6) 
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "environment") 
})

## Test 1.1.2 Misspecified Qform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when Qform is misspecified", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = Qform.mis, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3256744, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3278714, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.2859616, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Ygstar), abs(estimates["gcomp", ] - psi0.Ygstar))
})

## Test 1.1.3 Misspecified gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when gform is misspecified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3297466, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3297466, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.3293320, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Ygstar), abs(estimates["iptw", ] - psi0.Ygstar))
})

## Test 1.1.4 Misspecified gform (+ tmle.intercept + equal.mass + 20 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when gform is misspecified (20 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY), nbins = 20)
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3317774, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.2908742, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.3293320, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Ygstar), abs(estimates["iptw", ] - psi0.Ygstar))
})

######################################### 
## Test 1.2 h2o
#########################################
## Test 1.2.1 h2o.glm.wrapper (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with h2o, using h2o.glm.wrapper", {
  tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = nrow(dat_iidcontABinY),
                  h2olearner = c("h2o.glm.wrapper"), nbins = 10)
  require(h2oEnsemble)
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3268742, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3247688, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3291617, tolerance = 1e-6) 
  h2o::h2o.shutdown(prompt = FALSE)
})

## Test 1.2.2 h2o.glm.wrapper & h2o.randomForest.wrapper (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with h2o, using h2o.glm.wrapper", {
  tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = nrow(dat_iidcontABinY),
                  h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), nbins = 10)
  require(h2oEnsemble)
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3270933, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3154016, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3291424, tolerance = 1e-6) 
  h2o::h2o.shutdown(prompt = FALSE)
})

######################################### 
## Test 1.3 SuperLearner
#########################################
## Test 1.3.1 SL.glm, SL.step, SL.glm.interaction (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Binary Y) for continuous A with SuperLearner, using SL.glm, SL.step, SL.glm.interaction", {
  tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(dat_iidcontABinY),
                  g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"), nbins = 10)
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], 0.3279654, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3266310, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3266310, tolerance = 1e-6) 
})


# ---------------------------------------------------------------------------------
# Test 2. The TMLE/ IPTW/ GCOMP estimator that fits continuous density when f_gstar1 or f_gstar2 = NULL
# ---------------------------------------------------------------------------------
######################################### 
## Test 2.1 f_gstar1 = NULL & f_gstar2 = NULL
#########################################
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when f_gstar1 = NULL & f_gstar2 = NULL", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
  expect_message(tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                                              f_gstar1 = NULL, f_gstar2 = NULL, Qform = Qform.corr, savetime.fit.hbars = TRUE),
                 regexp = "skip g0 & gstar fitting procedure and directly set h_gstar_h_gN = 1 for each observation") 
  estimates <- tmleCom_res$EY_gstar1$estimates
  expect_equivalent(as.vector(estimates), rep(mean(dat_iidcontABinY$Y), 3))
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "NULL")
})

######################################### 
## Test 1.2 f_gstar1 = NULL & f_gstar2 != NULL
#########################################
test_that("fit TMLE estimator (Binary Y) for continuous A with speedglm, when f_gstar1 = NULL & f_gstar2 != NULL", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = NULL, f_gstar2 = f.gstar, Qform = Qform.corr, savetime.fit.hbars = TRUE)
  estimates1 <- tmleCom_res$EY_gstar1$estimates
  estimates2 <- tmleCom_res$EY_gstar2$estimates  # psi.o = 0.316274
  expect_equivalent(as.vector(estimates1), rep(mean(dat_iidcontABinY$Y), 3))
  expect_equal(estimates2["tmle", ], 0.3280404, tolerance = 0.001)   
  expect_equal(estimates2["iptw", ], 0.3278714, tolerance = 0.001)  
  expect_equal(estimates2["gcomp", ], 0.3293320, tolerance = 0.001)  
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "NULL")
  expect_type(tmleCom_res$EY_gstar2$h.g0_GenericModel, "environment")
})


# ---------------------------------------------------------------------------------
# Test 2. The TMLE/ IPTW/ GCOMP estimator for continuous Y 
# --------------------------------------------------------------------------------- 
data(sampleDat_iidcontAContY)
dat_iidcontAContY <- sampleDat_iidcontAContY$dat_iidcontABinY
head(dat_iidcontAContY)
psi0.Y <- sampleDat_iidcontAContY$psi0.Y  # 0.291398
psi0.Ygstar <- sampleDat_iidcontAContY$psi0.Ygstar  # 0.316274

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    set.seed(rndseed)
    A.mu <- 0.86 * data[,"W1"] + 0.41 * data[,"W2"] - 0.34 * data[,"W3"] + 0.93 * data[,"W4"]
    untrunc.A <- rnorm(n = nrow(data), mean = A.mu + shift.const, sd = 1)
    r.new.A <- exp(0.5 * shift.const * (untrunc.A - A.mu - shift.const / 2))
    trunc.A <- ifelse(r.new.A > trunc.const, untrunc.A - shift.const, untrunc.A)
    return(trunc.A)
  }
  return(f.gstar)
}
f.gstar <- define_f.gstar(shift = sampleDat_iidcontAContY$shift.val, truncBD = sampleDat_iidcontAContY$truncBD, 
                          rndseed = sampleDat_iidcontAContY$rndseed)

Qform.corr <- "Y ~ W1 + W2 + W3 + W4 + A" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g
gform.mis <- "A ~ W2 + W3"  # correct g

######################################### 
## Test 2.1 speed.glm & glm 
######################################### 
## Test 2.1.1 Correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (Continuous Y) for continuous A with speedglm, when Qform & gform are correctly specified", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontAContY))
  tmleCom_res <- tmleCommunity(data = dat_iidcontAContY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.316274 
  expect_equal(estimates["tmle", ], , tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], , tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], , tolerance = 1e-6) 
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "environment") 
})
