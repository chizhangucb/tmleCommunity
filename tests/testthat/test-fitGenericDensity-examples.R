# ---------------------------------------------------------------------------------
# TEST SET 4. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
# FOCUS ON using function fitGenericDensity in hbarDensityModel.R
# Fitting continuous exposure by  binning, conditional on covariates
# ---------------------------------------------------------------------------------
library(testthat)
library(h2o)
gvars$verbose <- TRUE

# ---------------------------------------------------------------------------------
# Test 1. The IPTW estimator that fits onr continuous A density
# ---------------------------------------------------------------------------------
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontABinY)
dat_iidcontABinY <- sampleDat_iidcontABinY$dat_iidcontABinY
head(dat_iidcontABinY)
psi0.Y <- sampleDat_iidcontABinY$psi0.Y  # 0.291398
psi0.Ygstar <- sampleDat_iidcontABinY$psi0.Ygstar  # 0.316274
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL)

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

test.fitGeneric.density <- function(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                                    h2olearner = "h2o.glm.wrapper", h2ometalearner = "h2o.glm.wrapper", 
                                    g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"), data,
                                    killh2o.atlast = TRUE) {
  tmleCom_Options(Qestimator = Qestimator, gestimator = gestimator, bin.method = "equal.mass", maxNperBin = nrow(data),
                  h2olearner = h2olearner, h2ometalearner = h2ometalearner, g.SL.library = g.SL.library)
  
  # -------------------------------------------------------------------------------------------
  # estimating h_g0 and h_gstar without bounding
  # -------------------------------------------------------------------------------------------
  h_gN <- fitGenericDensity(data, Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), f_gstar = NULL, lbound = 0)$h_gstar
  h_gstar <- fitGenericDensity(data, Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), f_gstar = f.gstar, lbound = 0)$h_gstar
  
  # -------------------------------------------------------------------------------------------
  max(h_gN); min(h_gN); print(max(h_gN)); print(min(h_gN))
  max(h_gstar); min(h_gstar); print(max(h_gstar)); print(min(h_gstar)) 
  
  trim_wt <- 1/0.005
  wts <- h_gstar / h_gN
  wts[is.nan(wts)] <- 0
  wts[wts > trim_wt] <- trim_wt
  summary(h_gstar/h_gN); print(summary(h_gstar/h_gN))
  summary(wts); print(summary(wts))
  
  iptw_untrimmed <- mean(data[,"Y"] * (h_gstar/h_gN))
  iptw_trimmed <- mean(data[,"Y"] * (wts))
  psi0 <- mean(data$Y.gstar)
  
  print("true psi0: " %+% psi0)
  print("iptw (untrimmed): " %+% round(iptw_untrimmed, 6))
  print("iptw (wts trimmed by " %+% trim_wt %+% "): " %+% round(iptw_trimmed, 6))
  
  if (any(c(getopt("Qestimator"), getopt("gestimator")) %in% "h2o__ensemble")) {
    if (killh2o.atlast) { h2o.shutdown(prompt = FALSE) }
  }
  # test 1:
  # checkTrue(abs(psi0 - 0.24023) < 10^-6)
  # test 2:
  # checkTrue(abs(iptw_untrimmed - 0.243353) < 10^-6) 
  # test 3:
  # checkTrue(abs(iptw_trimmed - 0.243353) < 10^-6)
  return(list(psi0 = psi0, iptw_untrimmed = iptw_untrimmed, iptw_trimmed = iptw_trimmed))
}

############################# 
## Test 1.1 speedglm__glm
#############################
# test.fitGeneric.density(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", data = Odata)
test_that("fit iptw estimator for continuous A with speedglm", {  # psi0 = 0.3196
  iptw_res <- test.fitGeneric.density(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", data = dat_iidcontABinY)
  expect_equal(iptw_res$iptw_untrimmed, 0.327871, tolerance = 0.001)
  expect_equal(iptw_res$iptw_trimmed, 0.327871, tolerance = 0.001)
})
# ------------------------------------------------------
# Benchmark for using speed.glm
# ------------------------------------------------------
# >   # summeas.gstar$getPsAsW.models()[[1]]$intrvls.width
# >   max(h_gN); min(h_gN)
# [1] 0.5227538
# [1] 0.0008420795
# >   max(h_gstar); min(h_gstar);
# [1] 0.4969672
# [1] 6.514282e-05
# > 
# >   wts <- h_gstar / h_gN
# >   wts[is.nan(wts)] <- 0
# >   # wts[wts > 200] <- 200
# > 
# >   summary(h_gstar/h_gN)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.03416   0.29970   0.57810   1.00300   1.17800   17.57000 
# >   summary(wts)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.03416   0.29970   0.57810   1.00300   1.17800   17.57000 
# >   (iptw <- mean(datO[,"Y"] * (wts)))
# [1] "true psi0: 0.3196"
# [1] "iptw (untrimmed): 0.327871"
# [1] "iptw (wts trimmed by 200): 0.327871"

############################# 
## Test 1.2 h2o
#############################
# test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", h2olearner = "h2o.glm.wrapper", 
#                         data = Odata, killh2o.atlast = TRUE)
test_that("fit iptw estimator for continuous A with only h2o.glm.wrapper algorithm in h2oEnsemble", {  # dat_iidcontABinY
  iptw_res <- test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", 
                                      h2olearner = "h2o.glm.wrapper", data = dat_iidcontABinY, killh2o.atlast = TRUE)
  require(h2oEnsemble)
  expect_equal(iptw_res$iptw_untrimmed, 0.324769, tolerance = 0.001)
  expect_equal(iptw_res$iptw_trimmed, 0.324769, tolerance = 0.001)
})
# ------------------------------------------------------
# Benchmark for h2olearner = "h2o.glm.wrapper" 
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.6031307
# [1] 0.0005458943
# >   max(h_gstar); min(h_gstar);  
# [1] 0.5520781  
# [1] 0.002283248  
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.1119     0.3038   0.5821    1.0010    1.1380    15.5500 
# >   summary(wts)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.1119     0.3038   0.5821    1.0010    1.1380    15.5500  
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 0.3196"  
# [1] "iptw (untrimmed): 0.324769"  
# [1] "iptw (wts trimmed by 200): 0.324769"  

# test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", 
#                         h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), data = Odata, killh2o.atlast = TRUE)
test_that("fit iptw estimator for continuous A with h2o.glm.wrapper & h2o.randomForest.wrapper algorithms in h2oEnsemble", {
  iptw_res <- test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", data = dat_iidcontABinY, 
                                      h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), killh2o.atlast = TRUE)
  require(h2oEnsemble)
  expect_equal(iptw_res$iptw_untrimmed, 0.31715, tolerance = 0.01)
  expect_equal(iptw_res$iptw_trimmed, 0.31715, tolerance = 0.01)
})
# ------------------------------------------------------
# Benchmark for h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper")  
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.6341132
# [1] 0.0006458569
# >   max(h_gstar); min(h_gstar);  
# [1] 0.7649564  
# [1] 0.00223  
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.09126   0.30790   0.59980   0.97690   1.11100   12.98000   
# >   summary(wts)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.09126   0.30790   0.59980   0.97690   1.11100   12.98000   
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 0.3196"  
# [1] "iptw (untrimmed): 0.31715"  
# [1] "iptw (wts trimmed by 200): 0.31715" 

# test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", data = Odata, 
#                         h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper", "h2o.gbm.wrapper"), killh2o.atlast = TRUE)
test_that("fit iptw estimator for continuous A with h2o.glm.wrapper, h2o.randomForest.wrapper & h2o.gbm.wrapper algorithms in h2oEnsemble", {
  iptw_res <- test.fitGeneric.density(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", data = dat_iidcontABinY,  
                                      killh2o.atlast = T, h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper", "h2o.gbm.wrapper"))
  require(h2oEnsemble)
  expect_equal(iptw_res$iptw_untrimmed, 0.322908, tolerance = 0.01)
  expect_equal(iptw_res$iptw_trimmed, 0.322908, tolerance = 0.01)
})
# ------------------------------------------------------
# Benchmark for h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper",  "h2o.gbm.wrapper")  
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.8315298
# [1] 0.0006624244
# >   max(h_gstar); min(h_gstar);  
# [1] 0.7282337
# [1] 0.002223962
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.08709   0.31560   0.60590   0.99530   1.13300   14.17000   
# >   summary(wts)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.08709   0.31560   0.60590   0.99530   1.13300   14.17000
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 0.3196"  
# [1] "iptw (untrimmed): 0.322908"  
# [1] "iptw (wts trimmed by 200): 0.322908"  

############################# 
## Test 1.3 SuperLearner
#############################
# test.fitGeneric.density(Qestimator = "SuperLearner", gestimator = "SuperLearner",  
#                         g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"), data = Odata)
test_that("fit iptw estimator for continuous A with SL.glm, SL.step & SL.glm.interaction algorithms in SuperLearner", {
  iptw_res <- test.fitGeneric.density(Qestimator = "SuperLearner", gestimator = "SuperLearner", data = dat_iidcontABinY,
                                      g.SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"))
  expect_equal(iptw_res$iptw_untrimmed, 0.3264449, tolerance = 0.005)
  expect_equal(iptw_res$iptw_trimmed, 0.3264449, tolerance = 0.005)
})
# ------------------------------------------------------
# Benchmark for using SuperLearner
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.5401792
# [1] 0.001040254
# >   max(h_gstar); min(h_gstar);  
# [1] 0.5112411  
# [1] 3.609393e-05  
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.01291   0.28630   0.59170   1.00100   1.18300   16.03000 
# >   summary(wts)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.01291   0.28630   0.59170   1.00100   1.18300   16.03000 
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 0.3196"  
# [1] "iptw (untrimmed): 0.3264449"  
# [1] "iptw (wts trimmed by 200): 0.3264449"  
