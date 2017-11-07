#***************************************************************************************
# Example 1: Hierarchical example, with one binary A and bianry Y 
# True ATE of the community-based treatment is approximately 0.103716
data(comSample.wmT.bA.bY_list)  # load the sample data 
comSample.wmT.bA.bY <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
N <- NROW(comSample.wmT.bA.bY)
Qform.corr <- "Y ~ E1 + E2 + W2 + W3 + A" # correct Q form
gform.corr <- "A ~ E1 + E2 + W1"  # correct g
#***************************************************************************************

#***************************************************************************************
# 1.1 Estimating the additive treatment effect (ATE) for two deterministic interventions
# (f_gstar1 = 1 vs f_gstar2 = 0) via community-level / individual-level analysis.
# speed.glm using correctly specified Qform, hform.g0 and hform.gstar;
#***************************************************************************************
# Setting global options that may be used in tmleCommunity(), e.g., using speed.glm
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)

# Community-level analysis without a pooled individual-level regression on outcome
tmleCom_wmT.bA.bY.1a_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = FALSE, 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)

# Examples of estimates under f_gstar1 = 1:
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar1$estimates
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar1$vars
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar1$CIs

# Examples of estimates under f_gstar0 = 0:
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar2$estimates
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar2$vars
tmleCom_wmT.bA.bY.1a_sglm$EY_gstar2$CIs

# Examples of estimates for ATE under f_gstar1 - f_gstar0:
tmleCom_wmT.bA.bY.1a_sglm$ATE$estimates
tmleCom_wmT.bA.bY.1a_sglm$ATE$vars
tmleCom_wmT.bA.bY.1a_sglm$ATE$CIs
head(tmleCom_wmT.bA.bY.1a_sglm$ATE$IC)

# Community-level analysis with a pooled individual-level regression on outcome
tmleCom_wmT.bA.bY.1b_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = TRUE, 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.1b_sglm$ATE$estimates

# Individual-level analysis with both individual-level outcome and treatment mechanisms
tmleCom_wmT.bA.bY.2_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "individual_level", communityID = "id", 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.2_sglm$ATE$estimates

# Failing to provide communityID will automatically set community.step to "NoCommunity"
tmleCom_wmT.bA.bY.NoC_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "individual_level", communityID = NULL, 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.NoC_sglm$ATE$estimates

#***************************************************************************************
# 1.2 Same as above but for different Qestimator and gestimator through tmleCom_Options()
# via community-level analysis with a pooled individual-level regression on outcome.
# (See more details in examples in tmleCom_Options())
#***************************************************************************************
# SuperLearner for both outcome and treatment (clever covariate) regressions
# using all parent nodes (of Y and A) as regressors (respectively)
require("SuperLearner")
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", 
                maxNperBin = N, SL.library = c("SL.glm", "SL.step", "SL.bayesglm"))
require("arm")
tmleCom_wmT.bA.bY.2_SL <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = TRUE, 
                Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
tmleCom_wmT.bA.bY.2_SL$ATE$estimates

# SuperLearner for outcome regressions and glm treatment regressions
# using all regressors in the correctly specified Qform and 
# all regressors in the misspecified hform.g0 and hform.gstar
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "glm__glm", 
                maxNperBin = N, SL.library = c("SL.mean", "SL.stepAIC", "SL.bayesglm"))
tmleCom_wmT.bA.bY.2_SL.glm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = TRUE,
                Qform = NULL, hform.g0 = "A ~ W1", hform.gstar = "A ~ E1 + W2")
tmleCom_wmT.bA.bY.2_SL.glm$ATE$estimates

#***************************************************************************************
# 1.3 Evaluating mean population outcome under static intervention A = 0
# with different community-level and individual-level weight choices 
#***************************************************************************************
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)
# weigh individuals in data equally & weigh community by its number of individuals
tmleCom_wmT.bA.bY.1a_w1 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                obs.wts = "equal.within.pop", community.wts = "size.community", 
                community.step = "community_level", communityID = "id")
tmleCom_wmT.bA.bY.1a_w1$EY_gstar1$estimates

# same as above but weigh individuals within the same community equally
tmleCom_wmT.bA.bY.1a_w2 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                obs.wts = "equal.within.community", community.wts = "size.community", 
                community.step = "community_level", communityID = "id")
tmleCom_wmT.bA.bY.1a_w2$EY_gstar1$estimates

# weigh individuals within the same community equally & weigh community equally
tmleCom_wmT.bA.bY.1a_w3 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                obs.wts = "equal.within.community", community.wts = "equal.community", 
                community.step = "community_level", communityID = "id")
tmleCom_wmT.bA.bY.1a_w3$EY_gstar1$estimates

#***************************************************************************************
# 1.4 Specifying user-supplied stochastic or deterministic intervention function
#***************************************************************************************
# Intervention function that will sample A with probability P(A=1) = prob.val
define_f.gstar <- function(prob.val, rndseed = NULL) {
  eval(prob.val) 
  f.gstar <- function(data, ...) {
    print("probability of selection: " %+% prob.val)
    rbinom(n = NROW(data), size = 1, prob = prob.val)
  }
  return(f.gstar)
}
# Stochastically set 50% of the population to A=1 
f.gstar_stoch.0.5 <- define_f.gstar(prob.val = 0.5)
# Deterministically set 100% of the population to A=1 
f.gstar_determ.1 <- define_f.gstar(prob.val = 1)
# Deterministically set 100% of the population to A=0
f.gstar_determ.0 <- define_f.gstar(prob.val = 0)

#***************************************************************************************
# 1.5 Equivalent ways of specifying user-supplied (static) intervention (f_gstar1 = 1)
#***************************************************************************************
# Alternative 1: via intervention functoin that sets every invidual's A to constant 1
tmleCom_wmT.bA.bY.1a_fgtar1 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar_determ.1, 
                community.step = "community_level", communityID = "id")
tmleCom_wmT.bA.bY.1a_fgtar1$EY_gstar1$estimates

# Alternative 2: by simply setting f_gstar1 to 1
tmleCom_wmT.bA.bY.1a_fgtar1 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L,
                community.step = "community_level", communityID = "id")

# Alternative 3: by setting f_gstar1 to a vector of 1's of length NROW(data)
tmleCom_wmT.bA.bY.1a_fgtar1 <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                f_gstar1 = rep(1L, NROW(comSample.wmT.bA.bY)), 
                community.step = "community_level", communityID = "id")

#***************************************************************************************
# 1.6 Running exactly the same estimator as 1.1 but using h_gstar/h_gN as a coviariate 
# in the targeting step (default to use weighted intercept-based TMLE)
#***************************************************************************************
# unweighted covariate-based TMLE
tmleCom_wmT.bA.bY.1a_covTMlE <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = FALSE, 
                TMLE.targetStep = "tmle.covariate",  # default as "tmle.intercept"
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.1a_covTMlE$ATE$estimates

#***************************************************************************************
# 1.7 Equivalent ways of specifying the regression formulae 
# (if Ynode is specified as "Y" and WEnodes = c("E1", "E2", "W1", "W2", "W3"))
#***************************************************************************************
# For outcome regression, the left side of Qform will be ignored if Ynode is specified,
# with dependent variable being set to Ynode.
Qform1 <- "Y ~ E1 + E2 + W2 + W3 + A" 
Qform2 <- "AnythingIsFine ~ E1 + E2 + W2 + W3 + A" 
Qform3 <- NULL  # since all parent nodes of Y will be used as regressors

# For treatment regressions, if hform.gstar unspecified, it uses the same regression
# formula as hform.g0 does. 
# Alternative 1: specify hform.g0 and hform.gstar respectively
hform.g0 <- "A ~ E1 + E2 + W1"
hform.gstar <- "A ~ E1 + E2 + W1"

# Alternative 2: specify hform.g0 only
hform.g0 <- "A ~ E1 + E2 + W1"
hform.gstar <- NULL

#***************************************************************************************
# Example 2: Non-hierarchical example, with one continuous A and continuous Y 
# True mean population outcome under stochastic intervention (specified below)
# is approximately 3.50856
data(indSample.iid.cA.cY_list)  # load the sample data 
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
true.shift <- indSample.iid.cA.cY_list$shift.val  # 2
true.truncBD <- indSample.iid.cA.cY_list$truncBD  # 10
N <- NROW(indSample.iid.cA.cY)
Qform.corr <- "Y ~ W1 + W2 + W3 + W4 + A" # correct Q form
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g
#***************************************************************************************

#***************************************************************************************
# 2.1 Specifying stochastic intervention function that could represent the true 
# shifted version of the current treatment mechanism
#***************************************************************************************
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
# correctly specified stochastic intervention with true shift value and truncated bound
f.gstar.corr <- define_f.gstar(shift = true.shift, truncBD = true.truncBD)
# Misspecified specified stochastic intervention 
f.gstar.mis <- define_f.gstar(shift = 5, truncBD = 8)

#***************************************************************************************
# 2.2 Estimating mean population outcome under different stochastic interventions
# speed.glm using correctly specified Qform, hform.g0 and hform.gstar
#***************************************************************************************
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)
# correctly specified stochastic intervention 
tmleind_iid.cA.cY_true.fgstar <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleind_iid.cA.cY_true.fgstar$EY_gstar1$estimates

# misspecified stochastic intervention
tmleind_iid.cA.cY_mis.fgstar <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.mis,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleind_iid.cA.cY_mis.fgstar$EY_gstar1$estimates

#***************************************************************************************
# 2.3 Same as above but using larger number of Monte-Carlo simulations
# using all parent nodes (of Y and A) as regressors (respectively)
#***************************************************************************************
# A will be sampled 10 times (for a total sample size of NROW(data)*10 under f_gstar1)
tmleind_iid.cA.cY_10MC <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr, n_MCsims = 10)
tmleind_iid.cA.cY_10MC$EY_gstar1$estimates

# A will be sampled 1000 times (for a total sample size of NROW(data)*1000 under f_gstar1)
tmleind_iid.cA.cY_1000MC <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr, n_MCsims = 1000)
tmleind_iid.cA.cY_1000MC$EY_gstar1$estimates

#***************************************************************************************
# 2.4 Same as above but printing out status messages 
#***************************************************************************************
tmleind_iid.cA.cY_10MC <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", verbose = TRUE, 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr, n_MCsims = 10)

#***************************************************************************************
# 2.5 Running exactly the same estimator as 2.1 but defining different values of bin cutoffs 
#***************************************************************************************
# using equal-length method with 10 bins 
tmleCom_Options(bin.method = "equal.len", nbins = 10, maxNperBin = N)
tmleind_iid.cA.cY_len <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleind_iid.cA.cY_len$EY_gstar1$estimates

# using combination of equal-length and equal-mass method with 20 bins 
tmleCom_Options(bin.method = "dhist", nbins = 20, maxNperBin = N)
tmleind_iid.cA.cY_dhist <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = f.gstar.corr,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleind_iid.cA.cY_dhist$EY_gstar1$estimates

#***************************************************************************************
# 2.6 Estimating the additive treatment effect (ATE) for two stochastic interventions
#***************************************************************************************
# Intervention function that will shift A by constant rate (shift.rate)
# (A special case of stochastic intervention with constant shift)
define_f.gstar <- function(shift.rate) {
  eval(shift.rate) 
  f.gstar <- function(data, ...) {
    print("rate of shift: " %+% shift.rate)
    data[, "A"] * shift.rate
  }
  return(f.gstar)
}
f.gstar_shift0.8 <- define_f.gstar(shift.rate = 0.8)
f.gstar_shift0.5 <- define_f.gstar(shift.rate = 0.6)

tmleCom_Options(bin.method = "equal.mass", nbins = 5, maxNperBin = N)
tmleind_iid.cA.cY_ATE <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"),
                f_gstar1 = f.gstar_shift0.8, f_gstar2 = f.gstar_shift0.5,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)

# ATE estimates for f_gstar1-f_gstar2:
tmleind_iid.cA.cY_ATE$ATE$estimates
tmleind_iid.cA.cY_ATE$ATE$vars
tmleind_iid.cA.cY_ATE$ATE$CIs   

#***************************************************************************************
# 2.7 Same as above but using a vector of user-supplied observation weights  
#***************************************************************************************
set.seed(12345)
obs.wts <- runif(n = NROW(indSample.iid.cA.cY), min = 0, max = 1)
tmleind_iid.cA.cY_own.obsWT <- 
  tmleCommunity(data = indSample.iid.cA.cY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("W1", "W2", "W3", "W4"), obs.wts = obs.wts,
                f_gstar1 = f.gstar_shift0.8, f_gstar2 = f.gstar_shift0.5,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleind_iid.cA.cY_own.obsWT$ATE$estimates
