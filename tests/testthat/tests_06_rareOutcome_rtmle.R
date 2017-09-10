# ---------------------------------------------------------------------------------
# TEST SET 5. TESTS FOR FITTING BINARY EXPOSURE A IN IID DATA WHERE OUTCOME IS RARE
# ---------------------------------------------------------------------------------
# AKA CASE-CONTROL STUDY

library(testthat)
gvars$verbose <- TRUE

######################################### 
## Test 1.1 Sample set with J = 1 (i.e., nCase/nControl = 1)
######################################### 
data(rareSamples_iidBinAY.J1)
# load("rareSamples_iidBinAY.J1.Rda")
dat_iidBinAY.rareJ1 <- rareSamples_iidBinAY.J1$dat_iidBinAY.rareJ1
obs.wt.J1 <- rareSamples_iidBinAY.J1$obs.wt.J1
q0 <- rareSamples_iidBinAY.J1$q0  # 0.013579
psi0.Y <- rareSamples_iidBinAY.J1$psi0.Y # 0.012662

Qform.corr <- "Y ~ W1 + W2*A + W3 + W4" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g
gform.mis <- "A ~ W2"  # incorrect g

## Test 1.1.1 Correct weights + correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for rare binary Y with correct weights (J=1), when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAY.rareJ1))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAY.rareJ1, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                               obs.wts = obs.wt.J1)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.0124744, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.0127704, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.0124322, tolerance = 1e-6) 
})

## Test 1.1.2 Incorrect weights + correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for rare binary Y with incorrect weights (J=1), when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAY.rareJ1))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAY.rareJ1, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                               obs.wts = NULL)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.2466571, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.2688000, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.2244907, tolerance = 1e-6) 
})

## Test 1.1.3 Correct weights + Misspecified Qform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for rare binary Y with correct weights (J=1), when Qform is misspecified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAY.rareJ1))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAY.rareJ1, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.mis, hform.g0 = gform.corr, hform.gstar = gform.corr,
                               obs.wts = obs.wt.J1)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.012548105, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.012770399, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.009324937, tolerance = 1e-6) 
})

######################################### 
## Test 1.2 Sample set with J = 2 (i.e., nCase/nControl = 1/2)
######################################### 
data(rareSamples_iidBinAY.J2)
# load("rareSamples_iidBinAY.J2.Rda")
dat_iidBinAY.rareJ2 <- rareSamples_iidBinAY.J2$dat_iidBinAY.rareJ2
obs.wt.J2 <- rareSamples_iidBinAY.J2$obs.wt.J2
q0 <- rareSamples_iidBinAY.J2$q0  # 0.013579
psi0.Y <- rareSamples_iidBinAY.J2$psi0.Y # 0.012662

Qform.corr <- "Y ~ W1 + W2*A + W3 + W4" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g
gform.mis <- "A ~ W2"  # incorrect g

## Test 1.1.1 Correct weights + correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for rare binary Y with correct weights (J=2), when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAY.rareJ2))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAY.rareJ2, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                               obs.wts = obs.wt.J2)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.01260533, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.01273315, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.01234571, tolerance = 1e-6) 
})

## Test 1.1.2 Incorrect weights + correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator for rare binary Y with correct weights (J=2), when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAY.rareJ2))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAY.rareJ2, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                               obs.wts = NULL)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.2091487, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.2197760, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.1944231, tolerance = 1e-6) 
})
