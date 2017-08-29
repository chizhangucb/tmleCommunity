# ---------------------------------------------------------------------------------
# TEST SET 1. TESTS FOR FITTING BINARY EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
# Fitting binary exposure by logistic regression
# ---------------------------------------------------------------------------------
library(testthat)
gvars$verbose <- TRUE

# ---------------------------------------------------------------------------------
# Test 1. The TMLE/ IPTW/ GCOMP estimator for binary Y 
# --------------------------------------------------------------------------------- 
# A is bernoulli for each observation being a function of (W1, W2, W3, W4);
data(sampleDat_iidBinABinY)
# load("sampleDat_iidBinABinY.Rda")
dat_iidBinABinY <- sampleDat_iidBinABinY$dat_iidBinABinY
dat_iidBinABinY <- dat_iidBinABinY[, c("W1", "W2", "W3", "W4", "A", "Y")]
head(dat_iidBinABinY)
psi0.Y <- sampleDat_iidBinABinY$psi0.Y  # 0.348242

Qform.corr <- "Y ~ W1 + W2*A + W3 + W4" # correct Q
Qform.mis <- "Y ~ W3 + A" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 * W4"  # correct g
gform.mis <- "A ~ W2 + W3"  # incorrect g

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
## Test 1.1.1 Correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with speedglm, when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3513866, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.3506032, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3507349, tolerance = 1e-6) 
  expect_type(tmleCom_res$EY_gstar1$h.g0_GenericModel, "environment") 
})

## Test 1.1.2 Misspecified Qform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with speedglm, when Qform is misspecified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.mis, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3510578, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3506032, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.4494514, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["gcomp", ] - psi0.Y))
})

## Test 1.1.3 Misspecified gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with speedglm, when gform is misspecified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinABinY))
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3501777, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.4129737, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.3507349, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["iptw", ] - psi0.Y))
})

## Test 1.1.4 Misspecified gform (+ tmle.intercept + equal.mass + 20 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with speedglm, when gform is misspecified (20 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinABinY), nbins = 20)
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3501777, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.4129737, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.3507349, tolerance = 1e-6) 
  expect_lt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["iptw", ] - psi0.Y))
})

## Test 1.1.5  Misspecified Qform & gform with only main terms (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with speedglm, when both Qform and gform only use main terms", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinABinY), nbins = 10)
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = "Y ~ A + W1 + W2 + W3 + W4", 
                               hform.g0 = "A ~ W1 + W2 + W3 + W4", hform.gstar = "A ~ W1 + W2 + W3 + W4")
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3508978, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3508871, tolerance = 1e-6) 
  expect_equal(estimates["gcomp", ], 0.3470415, tolerance = 1e-6) 
  expect_gt(abs(estimates["tmle", ] - psi0.Y), abs(estimates["gcomp", ] - psi0.Y))
})

######################################### 
## Test 1.2 h2o
#########################################
## Test 1.2.1 h2o.glm.wrapper with only main terms (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for binary A with h2o, using h2o.glm.wrapper (10 bins),
           when both Qform and gform only use main terms", {
  tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = nrow(dat_iidBinABinY),
                  h2olearner = c("h2o.glm.wrapper"), nbins = 10)
  require(h2oEnsemble)           
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3516327, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3541630, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3166018, tolerance = 1e-6) 
  h2o.shutdown(prompt = FALSE)
  # If Qform is NOT SPECIFIED, all available main terms in (A, W, E) will be included into predictor set
  # If hform.g0 is NOT SPECIFIED, all available main terms in (W, E) will be included into predictor set, so is hform.gstar
})

######################################### 
## Test 1.3 SuperLearner
#########################################
## Test 1.3.1 SL.glm, SL.step, SL.glm.interaction with only main terms (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (binary Y) for continuous A with SuperLearner, using SL.glm, SL.step, SL.glm.interaction (10 bins),
           when both Qform and gform only use main terms", {
  tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(dat_iidBinABinY),
                  g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"), nbins = 10)
  tmleCom_res <- tmleCommunity(data = dat_iidBinABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.348242 
  expect_equal(estimates["tmle", ], 0.3512338, tolerance = 1e-6) 
  expect_equal(estimates["iptw", ], 0.3507137, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.3490567, tolerance = 1e-6) 
})


# ---------------------------------------------------------------------------------
# Test 2. The TMLE/ IPTW/ GCOMP estimator for continuous Y 
# --------------------------------------------------------------------------------- 
data(sampleDat_iidBinAContY)
# load("sampleDat_iidBinAContY.Rda")
dat_iidBinAContY <- sampleDat_iidBinAContY$dat_iidBinAContY
dat_iidBinAContY <- dat_iidBinAContY[, c("W1", "W2", "W3", "W4", "A", "Y")]
head(dat_iidBinAContY)
psi0.Y <- sampleDat_iidBinAContY$psi0.Y  # 2.8002673

Qform.corr <- "Y ~ A + W1 + W2 + W3 + W4" # correct Q
Qform.mis <- "Y ~ A + W2" # incorrect Q
gform.corr <- "A ~ W1 + W2 + W3 * W4"  # correct g
gform.mis <- "A ~ W2 + W3"  # incorrect g

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
## Test 1.1.1 Correctly specified Qform & gform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (continuous Y) for binary A with speedglm, when Qform & gform are correctly specified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAContY))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAContY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.8002673 
  expect_equal(estimates["tmle", ], 2.7893912, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 2.7760817, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 2.7669255, tolerance = 1e-6) 
})

## Test 1.1.2  Misspecified Qform (+ tmle.intercept + equal.mass + 10 nbins)
test_that("fit TMLE estimator (continuous Y) for binary A with speedglm, when Qform is misspecified (10 bins)", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidBinAContY))
  tmleCom_res <- tmleCommunity(data = dat_iidBinAContY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                               f_gstar1 = 1, f_gstar2 = 0, Qform = Qform.mis, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 2.8002673 
  expect_equal(estimates["tmle", ], 2.7941109, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 2.7760817, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 3.5052616, tolerance = 1e-6) 
})
