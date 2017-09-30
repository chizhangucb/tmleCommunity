# ---------------------------------------------------------------------------------
# TEST SET 5. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN HIERARCHICAL DATA
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by logistic regression
# ---------------------------------------------------------------------------------
source("tmleCommunity/tests/dataGeneration/get.cluster.dat.Acont.R")
library(testthat)
library(Hmisc)
gvars$verbose <- TRUE

# ---------------------------------------------------------------------------------
# Test 1. The TMLE/ IPTW/ GCOMP estimator for binary Y 
# --------------------------------------------------------------------------------- 

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
truncBD <- 5
shift.val <- 1

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    set.seed(rndseed)
    A.mu <- - 1.2  + 0.8 * data[,"E1"] + 0.21 * data[,"E2"] + 3 * data[,"W1"] - 0.7 * data[,"W2"] + 0.3 * data[,"W3"]
    untrunc.A <- rnorm(n = nrow(data), mean = A.mu + shift.const, sd = 1)
    r.new.A <- exp(1.5 * shift.const * (untrunc.A - A.mu - shift.const / 4))
    trunc.A <- ifelse(r.new.A > trunc.const, untrunc.A - shift.const, untrunc.A)
    return(trunc.A)
  }
  return(f.gstar)
}
f.gstar <- define_f.gstar(shift = shift.val, truncBD = truncBD, rndseed = NULL)

Qform.corr <- "Y ~ E1 + E2 + W1 + W2 + W3 + A" # correct Q
Qform.mis <- "Y ~ E1 + W1" # incorrect Q
gform.corr <- "A ~ E1 + E2 + W1 + W2 + W3"  # correct g
gform.mis <- "A ~ E1 + W3"  # incorrect g

bootstrap.TMLE <- function(nBoot, J, n.ind, truncBD, shift.val, nbins = 5, Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                           f.gstar, lbound = 0.005, working.model = TRUE, Qform = NULL, gform = NULL, rndseed = NULL) {
  
  est_tmle1a <- est_tmle1b <- est_tmle2 <- est_tmlePer <- NULL
  var_tmle1a <- var_tmle1b <- var_tmle2 <- var_tmlePer <- NULL
  CI_tmle1a <- CI_tmle1b <- CI_tmle2 <- CI_tmlePer <- NULL
  pval_tmle1a <- pval_tmle1b <- pval_tmle2 <- pval_tmlePer <- NULL
  error_dat <- error_1a <- error_1b <- error_2 <- error_Per <- 0
  
  set.seed(rndseed)
  r <- 1
  while (r < nBoot + 1) {
    message("##################################################################")
    message("################### Doing the " %+% r %+% "th bootstrap ###################")
    message("##################################################################\n")
    
    # Generate a sample of hierarchical data, when working.model fails
    data <- try(get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = NULL, truncBD = truncBD, shift.val = shift.val, 
                                  is.Y.bin = TRUE, working.model = working.model, timevarying = TRUE, n.ind.fix = FALSE))
    if (inherits(data, "try-error")) { error_dat <- error_dat + 1; next } # skip iterations with errors
    
    tmleCom_Options(Qestimator = Qestimator, gestimator = gestimator, maxNperBin = nrow(data), nbins = nbins)
    
    # TMLE-Ia
    tmleCom_1a_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "community_level", working.model = FALSE, f_gstar1 = f.gstar,
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound))
    if (inherits(tmleCom_1a_size, "try-error")) { error_1a <- error_1a + 1; next } # skip iterations with errors
    est_tmle1a <- rbind(est_tmle1a, tmleCom_1a_size$EY_gstar1$estimates[, 1])
    var_tmle1a <- rbind(var_tmle1a, tmleCom_1a_size$EY_gstar1$vars[, 1])
    CI_tmle1a <- rbind(CI_tmle1a, c(tmleCom_1a_size$EY_gstar1$CIs[, 1], tmleCom_1a_size$EY_gstar1$CIs[, 2]))
    pval_tmle1a <- rbind(pval_tmle1a, tmleCom_1a_size$EY_gstar1$pval[, 1])
    
    # TMLE-Ib
    tmleCom_1b_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "individual_level", working.model = FALSE,  f_gstar1 = f.gstar, 
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound))
    if (inherits(tmleCom_1b_size, "try-error")) { error_1b <- error_1b + 1; next } # skip iterations with errors
    est_tmle1b <- rbind(est_tmle1b, tmleCom_1b_size$EY_gstar1$estimates[, 1])
    var_tmle1b <- rbind(var_tmle1b, tmleCom_1b_size$EY_gstar1$vars[, 1])
    CI_tmle1b <- rbind(CI_tmle1b, c(tmleCom_1b_size$EY_gstar1$CIs[, 1], tmleCom_1b_size$EY_gstar1$CIs[, 2]))
    pval_tmle1b <- rbind(pval_tmle1b, tmleCom_1b_size$EY_gstar1$pval[, 1])
    
    # TMLE-II
    tmleCom_ii_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "individual_level", working.model = TRUE, f_gstar1 = f.gstar,
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound))
    if (inherits(tmleCom_ii_size, "try-error")) { error_2 <- error_2 + 1; next } # skip iterations with errors
    est_tmle2 <- rbind(est_tmle2, tmleCom_ii_size$EY_gstar1$estimates[, 1])
    var_tmle2 <- rbind(var_tmle2, tmleCom_ii_size$EY_gstar1$vars[, 1])
    CI_tmle2 <- rbind(CI_tmle2, c(tmleCom_ii_size$EY_gstar1$CIs[, 1], tmleCom_ii_size$EY_gstar1$CIs[, 2]))
    pval_tmle2 <- rbind(pval_tmle2, tmleCom_ii_size$EY_gstar1$pval[, 1])
    
    # TMLE-Per
    tmleCom_per_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                          obs.wts = "equal.within.community", community.wts = "equal.community", communityID = "id",
                                          community.step = "perCommunity", working.model = TRUE, f_gstar1 = f.gstar,
                                          Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound))
    if (inherits(tmleCom_per_size, "try-error")) { error_Per <- error_Per + 1; next } # skip iterations with errors
    est_tmlePer <- rbind(est_tmlePer, tmleCom_per_size$EY_gstar1$estimates[, 1])
    var_tmlePer <- rbind(var_tmlePer, tmleCom_per_size$EY_gstar1$vars[, 1])
    CI_tmlePer <- rbind(CI_tmlePer, c(tmleCom_per_size$EY_gstar1$CIs[, 1], tmleCom_per_size$EY_gstar1$CIs[, 2]))
    pval_tmlePer <- rbind(pval_tmlePer, tmleCom_per_size$EY_gstar1$pval[, 1])
    
    # Move forward once all four TMLEs have been calculated in this iteration
    r <- r + 1
  }
  return(list(est_tmle1a = est_tmle1a, est_tmle1b = est_tmle1b, est_tmle2 = est_tmle2, est_tmlePer = est_tmlePer,
              var_tmle1a = var_tmle1a, var_tmle1b = var_tmle1b, var_tmle2 = var_tmle2, var_tmlePer = var_tmlePer,
              CI_tmle1a = CI_tmle1a, CI_tmle1b = CI_tmle1b, CI_tmle2 = CI_tmle2, CI_tmlePer = CI_tmlePer,
              pval_tmle1a = pval_tmle1a, pval_tmle1b = pval_tmle1b, pval_tmle2 = pval_tmle2, pval_tmlePer = pval_tmlePer,
              error = c(error_dat, error_1a, error_1b, error_2, error_Per)))
}

# Want to test combination of J = c(30, 300, 500, 1000) and n = c(30, 300, 500) whne conducting 1000 simulations

#######################################################
############# Whne working model holds ################
#######################################################
## Test 1.1.1 Correctly specified Qform & gform
## Test 1.1.1.1 The number of communities = 30 & the number of individual per community = 30
tmle.wmT.Qc.gc.J30.n30.T <- 
  system.time(tmle.wmT.Qc.gc.J30.n30 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 30, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                       Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                       lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J30.n30, tmle.wmT.Qc.gc.J30.n30.T, file = "tmle.wmT.Qc.gc.J30.n30.Rda")

## Test 1.1.1.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmT.Qc.gc.J300.n300.T <- 
  system.time(tmle.wmT.Qc.gc.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 200, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J300.n300, tmle.wmT.Qc.gc.J300.n300.T, file = "tmle.wmT.Qc.gc.J300.n300.Rda")

## Test 1.1.1.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmT.Qc.gc.J500.n300.T <-
  system.time(tmle.wmT.Qc.gc.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J500.n300, tmle.wmT.Qc.gc.J500.n300.T, file = "tmle.wmT.Qc.gc.J500.n300.Rda")

## Test 1.1.1.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmT.Qc.gc.J1000.n300.T <- 
  system.time(tmle.wmT.Qc.gc.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J1000.n300, tmle.wmT.Qc.gc.J1000.n300.T, file = "tmle.wmT.Qc.gc.J1000.n300.Rda")

## Test 1.1.1.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmT.Qc.gc.J300.n500.T <- 
  system.time(tmle.wmT.Qc.gc.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J300.n500, tmle.wmT.Qc.gc.J300.n500.T, file = "tmle.wmT.Qc.gc.J300.n500.Rda")

## Test 1.1.1.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmT.Qc.gc.J500.n500.T <-
  system.time(tmle.wmT.Qc.gc.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qc.gc.J500.n500, tmle.wmT.Qc.gc.J500.n500.T, file = "tmle.wmT.Qc.gc.J500.n500.Rda")

## Test 1.1.2 Correctly specified Qform & Misspecified gform
## Test 1.1.2.1 The number of communities = 30 & the number of individual per community = 300
tmle.wmT.Qc.gm.J30.n300.T <- 
  system.time(tmle.wmT.Qc.gm.J30.n300 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                        Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                        lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J30.n300, tmle.wmT.Qc.gm.J30.n300.T, file = "tmle.wmT.Qc.gm.J30.n300.Rda")

## Test 1.1.2.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmT.Qc.gm.J300.n300.T <-
  system.time(tmle.wmT.Qc.gm.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J300.n300,tmle.wmT.Qc.gm.J300.n300.T, file = "tmle.wmT.Qc.gm.J300.n300.Rda")

## Test 1.1.2.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmT.Qc.gm.J500.n300.T <- 
  system.time(tmle.wmT.Qc.gm.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J500.n300, tmle.wmT.Qc.gm.J500.n300.T, file = "tmle.wmT.Qc.gm.J500.n300.Rda")

## Test 1.1.2.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmT.Qc.gm.J1000.n300.T <- 
  system.time(tmle.wmT.Qc.gm.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J1000.n300, tmle.wmT.Qc.gm.J1000.n300.T, file = "tmle.wmT.Qc.gm.J1000.n300.Rda")

## Test 1.1.2.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmT.Qc.gm.J300.n500.T <- 
  system.time(tmle.wmT.Qc.gm.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J300.n500, tmle.wmT.Qc.gm.J300.n500.T, file = "tmle.wmT.Qc.gm.J300.n500.Rda")

## Test 1.1.2.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmT.Qc.gm.J500.n500.T <-
  system.time(tmle.wmT.Qc.gm.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmT.Qc.gm.J500.n500, tmle.wmT.Qc.gm.J500.n500.T, file = "tmle.wmT.Qc.gm.J500.n500.Rda")

## Test 1.1.3 Misspecified Qform & Correctly specified gform
## Test 1.1.3.1 The number of communities = 30 & the number of individual per community = 300
tmle.wmT.Qm.gc.J30.n300.T <- 
  system.time(tmle.wmT.Qm.gc.J30.n300 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                        Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                        lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J30.n300, tmle.wmT.Qm.gc.J30.n300.T, file = "tmle.wmT.Qm.gc.J30.n300.Rda")

## Test 1.1.3.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmT.Qm.gc.J300.n300.T <- 
  system.time(tmle.wmT.Qm.gc.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J300.n300, tmle.wmT.Qm.gc.J300.n300.T, file = "tmle.wmT.Qm.gc.J300.n300.Rda")

## Test 1.1.3.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmT.Qm.gc.J500.n300.T <- 
  system.time(tmle.wmT.Qm.gc.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J500.n300, tmle.wmT.Qm.gc.J500.n300.T, file = "tmle.wmT.Qm.gc.J500.n300.Rda")

## Test 1.1.3.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmT.Qm.gc.J1000.n300.T <- 
  system.time(tmle.wmT.Qm.gc.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J1000.n300, tmle.wmT.Qm.gc.J1000.n300.T, file = "tmle.wmT.Qm.gc.J1000.n300.Rda")

## Test 1.1.3.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmT.Qm.gc.J300.n500.T <- 
  system.time(tmle.wmT.Qm.gc.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J300.n500, tmle.wmT.Qm.gc.J300.n500.T, file = "tmle.wmT.Qm.gc.J300.n500.Rda")

## Test 1.1.3.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmT.Qm.gc.J500.n500.T <- 
  system.time(tmle.wmT.Qm.gc.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = TRUE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmT.Qm.gc.J500.n500, tmle.wmT.Qm.gc.J500.n500.T, file = "tmle.wmT.Qm.gc.J500.n500.Rda")
