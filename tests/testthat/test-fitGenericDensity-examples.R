context("Test fitGenericDensity")

# ---------------------------------------------------------------------------------
# Focus on testing the IPTW estimator that fits on continuous A density
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
# ---------------------------------------------------------------------------------
`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
N <- nrow(indSample.iid.cA.cY)
psi0.Y <- indSample.iid.cA.cY_list$psi0.Y  # 0.291398
psi0.Ygstar <- indSample.iid.cA.cY_list$psi0.Ygstar  # 0.316274
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))

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
f.gstar <- define_f.gstar(shift = indSample.iid.cA.cY_list$shift.val,
                          truncBD = indSample.iid.cA.cY_list$truncBD)

test.fitGeneric.density <- function(data, killh2o.atEnd = TRUE, estimator = "speedglm__glm", 
                                    h2olearner = "h2o.glm.wrapper", h2ometalearner = "h2o.glm.wrapper", 
                                    SL.library = c("SL.glm", "SL.step", "SL.glm.interaction")) {
  tmleCom_Options(gestimator = estimator, bin.method = "equal.mass", 
                  maxNperBin = NROW(data), nbins = 10, h2olearner = h2olearner, 
                  h2ometalearner = h2ometalearner, SL.library = SL.library)
  
  # -------------------------------------------------------------------------------------------
  # estimating h_g0 and h_gstar without bounding
  # -------------------------------------------------------------------------------------------
  h_gN <- fitGenericDensity(data, Anodes = nodes$Anodes, Wnodes = nodes$WEnodes, 
                            f_gstar = NULL, lbound = 0)$h_gstar
  h_gstar <- fitGenericDensity(data, Anodes = nodes$Anodes, Wnodes = nodes$WEnodes,
                               f_gstar = f.gstar, lbound = 0)$h_gstar
  
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
  
  if (tmleCommunity:::getopt("Qestimator") %in% "h2o__ensemble") {
    if (cleanh2o.atEnd) { h2o.removeAll() }
  }
  # test 1:
  # expect_true(abs(psi0 - 3.497218) < 10^-6)
  # test 2:
  # expect_true(abs(iptw_untrimmed - 3.444627) < 10^-6) 
  # test 3:
  # expect_true(abs(iptw_trimmed - 3.444627) < 10^-6)
  return(list(psi0 = psi0, iptw_untrimmed = iptw_untrimmed, iptw_trimmed = iptw_trimmed))
}

test_that("fit iptw estimator for continuous A with speedglm", {  # psi0 = 3.497218
  set.seed(12345)
  iptw_res <- test.fitGeneric.density(data = indSample.iid.cA.cY, estimator = "speedglm__glm")
  expect_equal(iptw_res$iptw_untrimmed, 3.445208, tolerance = 0.01)
  expect_equal(iptw_res$iptw_trimmed, 3.445208, tolerance = 0.01)
})
# ------------------------------------------------------
# Benchmark for using speed.glm (It changes each time)
# -----------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.5227538
# [1] 0.0008420795
# >   max(h_gstar); min(h_gstar);
# [1] 0.8379246
# [1] 7.289328e-07
# > 
# >   wts <- h_gstar / h_gN
# >   wts[is.nan(wts)] <- 0
# >   # wts[wts > 200] <- 200
# > 
# >   summary(h_gstar/h_gN)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.000292  0.038228  0.404475  0.996869  1.551241  5.441820
# >   summary(wts)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#   0.000292  0.038228  0.404475  0.996869  1.551241  5.441820
# >   (iptw <- mean(datO[,"Y"] * (wts)))
# [1] "true psi0: 3.497218"
# [1] "iptw (untrimmed): 3.445208"
# [1] "iptw (wts trimmed by 200): 3.445208"

test_that("fit iptw for cont A with SL.glm, SL.step & SL.glm.interaction in SuperLearner", {
  require("SuperLearner")
  set.seed(12345)
  iptw_res <- 
    test.fitGeneric.density(data = indSample.iid.cA.cY, estimator = "SuperLearner",
                            SL.library = c("SL.glm", "SL.step", "SL.glm.interaction"))
  expect_equal(iptw_res$iptw_untrimmed, 3.440705, tolerance = 0.05)
  expect_equal(iptw_res$iptw_trimmed, 3.440705, tolerance = 0.05)
})
# ------------------------------------------------------
# Benchmark for using SuperLearner (It changes each time)
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
#   0.000298  0.033324   0.413184   0.996053  1.567145   6.682385 
# >   summary(wts)  
#      Min.    1st Qu.    Median      Mean     3rd Qu.      Max.  
#   0.000298  0.033324   0.413184   0.996053  1.567145   6.682385 
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 3.497218"  
# [1] "iptw (untrimmed): 3.440705"  
# [1] "iptw (wts trimmed by 200): 3.440705"  

test_that("fit iptw estimator for continuous A with only h2o.glm.wrapper algorithm in h2oEnsemble", {
  require(h2oEnsemble)
  set.seed(12345)
  iptw_res <- test.fitGeneric.density(data = indSample.iid.cA.cY, estimator = "h2o__ensemble", 
                                      h2olearner = "h2o.glm.wrapper")
  expect_equal(iptw_res$iptw_untrimmed, 3.415603, tolerance = 0.05)
  expect_equal(iptw_res$iptw_trimmed, 3.415603, tolerance = 0.05)
})
# ------------------------------------------------------
# Benchmark for h2olearner = "h2o.glm.wrapper" (It changes each time)
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.6031307
# [1] 0.0005458943
# >   max(h_gstar); min(h_gstar);  
# [1] 1  
# [1] 0.0003175793  
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.     1st Qu.    Median      Mean      3rd Qu.     Max.  
#   0.005515  0.079717   0.386387    0.993738   1.373648   7.425610 
# >   summary(wts)  
#      Min.     1st Qu.    Median      Mean      3rd Qu.     Max.  
#   0.005515  0.079717   0.386387    0.993738   1.373648   7.425610 
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 3.497218"  
# [1] "iptw (untrimmed): 3.415603"  
# [1] "iptw (wts trimmed by 200): 3.415603"  

test_that("fit iptw estimator for continuous A with h2o.glm.wrapper & " %+% 
            "h2o.randomForest.wrapper algorithms in h2oEnsemble", {
  require(h2oEnsemble)           
  set.seed(12345)
  iptw_res <- test.fitGeneric.density(data = indSample.iid.cA.cY, estimator = "h2o__ensemble", 
                                      h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"))
  expect_equal(iptw_res$iptw_untrimmed, 3.378611, tolerance = 0.1)
  expect_equal(iptw_res$iptw_trimmed, 3.378611, tolerance = 0.1)
})

# ------------------------------------------------------
# Benchmark for h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper")  
# (It changes each time)
# ------------------------------------------------------
# >   max(h_gN); min(h_gN)
# [1] 0.667917
# [1] 0.0006049437
# >   max(h_gstar); min(h_gstar);  
# [1] 1  
# [1] 2.031023e-05  
# >   
# >   wts <- h_gstar / h_gN  
# >   wts[is.nan(wts)] <- 0  
# >   # wts[wts > 200] <- 200  
# >   
# >   summary(h_gstar/h_gN)  
#      Min.   1st Qu.    Median      Mean    3rd Qu.      Max.  
#   0.000166 0.102947  0.405738   0.982751  1.351620   8.286293 
# >   summary(wts)  
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
#   0.000166 0.102947  0.405738   0.982751  1.351620   8.286293
# >   (iptw <- mean(datO[,"Y"] * (wts)))  
# [1] "true psi0: 3.497218"  
# [1] "iptw (untrimmed): 3.378611"  
# [1] "iptw (wts trimmed by 200): 3.378611" 

#test_that("fit iptw estimator for continuous A with h2o.glm.wrapper, h2o.randomForest.wrapper & h2o.gbm.wrapper", {
#  require(h2oEnsemble)  
#  set.seed(12345)
#  iptw_res <- 
#    test.fitGeneric.density(data = indSample.iid.cA.cY, estimator = "h2o__ensemble", 
#                            h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper", "h2o.gbm.wrapper"))
#  expect_equal(iptw_res$iptw_untrimmed, **, tolerance = 0.01)
#  expect_equal(iptw_res$iptw_trimmed, **, tolerance = 0.01)
#})
