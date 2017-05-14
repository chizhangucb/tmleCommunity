# ---------------------------------------------------------------------------------
# TEST SET 4. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
# Show how to use FUNCTION fitGenericDensity for fitting arbitrary density
library(testthat)
gvars$verbose <- TRUE

# ---------------------------------------------------------------------------------
# Test 1. The IPTW estimator that fits one continuous A density (Same as tests in tests_04_fitIPTW_hbarDensityModel.R)
# ---------------------------------------------------------------------------------
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
data(sampleDat_iidcontABinY)
dat_iidcontABinY <- sampleDat_iidcontABinY$dat_iidcontABinY
psi0.Y <- sampleDat_iidcontABinY$psi0.Y  # 0.291398
psi0.Ygstar <- sampleDat_iidcontABinY$psi0.Ygstar  # 0.316274
nodes <- list(Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), Enodes = NULL, StratifyInd = NULL, YnodeDet = NULL)

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
f.gstar <- define_f.gstar(shift = sampleDat_iidcontABinY$shift.val, truncBD = sampleDat_iidcontABinY$truncBD, rndseed = sampleDat_iidcontABinY$rndseed)

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
  
  trim_wt <- 1/0.005
  wts <- h_gstar / h_gN
  wts[is.nan(wts)] <- 0
  wts[wts > trim_wt] <- trim_wt
  iptw_untrimmed <- mean(data[,"Y"] * (h_gstar/h_gN))
  iptw_trimmed <- mean(data[,"Y"] * (wts))
  psi0 <- mean(data$Y.gstar)
  
  if (any(c(getopt("Qestimator"), getopt("gestimator")) %in% "h2o__ensemble")) {
    if (killh2o.atlast) { h2o.shutdown(prompt = FALSE) }
  }
  return(list(psi0 = psi0, iptw_untrimmed = iptw_untrimmed, iptw_trimmed = iptw_trimmed))
}

test_that("fit iptw estimator for continuous A with speedglm", {  # psi0 = 0.3196
  iptw_res <- test.fitGeneric.density(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", data = dat_iidcontABinY)
  expect_equal(iptw_res$iptw_untrimmed, 0.327871, tolerance = 0.001)
  expect_equal(iptw_res$iptw_trimmed, 0.327871, tolerance = 0.001)
})


# ---------------------------------------------------------------------------------
# Test 2. 
# ---------------------------------------------------------------------------------

