#***************************************************************************************
# Example 1: Hierarchical example,  with one binary A and bianry Y 
#***************************************************************************************
# True ATE of the community-based treatment is approximately 0.103716
data(comSample.wmT.bA.bY_list)  # load the sample data 
comSample.wmT.bA.bY <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
Qform.corr <- "Y ~ E1 + E2 + W2 + W3 + A" # correct Q form
gform.corr <- "A ~ E1 + E2 + W1"  # correct g

### 1.1 Community-level analysis without a pooled individual-level regression on outcome
## 1.1.1 speed.glm using correctly specified Qform, hform.g0 and hform.gstar
# Setting global options that may be used in tmleCommunity(), e.g., using speed.glm
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                bin.method = "equal.mass", maxNperBin = nrow(data))

# Two weights choice "equal.within.pop" and "size.community"
tmleCom_wmT.bA.bY_Qcgc.1a <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), obs.wts = "equal.within.pop",
                community.step = "community_level", community.wts = "size.community", 
                communityID = "id", pooled.Q = FALSE, f_gstar1 = 1, f_gstar2 = 0,
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)

# examples of estimates under f_gstar1 = 1:
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar1$estimates
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar1$vars
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar1$CIs

# examples of estimates under f_gstar0 = 0:
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar2$estimates
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar2$vars
tmleCom_wmT.bA.bY_Qcgc.1a$EY_gstar2$CIs

# examples of estimates for ATE under f_gstar1 - f_gstar0:
tmleCom_wmT.bA.bY_Qcgc.1a$ATE$estimates
tmleCom_wmT.bA.bY_Qcgc.1a$ATE$vars
tmleCom_wmT.bA.bY_Qcgc.1a$ATE$CIs

## 1.2.1 SuperLearner using all parent nodes (of Y and A) as regressors (respectively)
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", 
                bin.method = "equal.mass", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.step", "SL.glm.interaction", "SL.bayesglm"))
tmleCom_wmT.bA.bY_SL.1a <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), obs.wts = "equal.within.pop",
                community.step = "community_level", community.wts = "size.community", 
                communityID = "id", pooled.Q = FALSE, f_gstar1 = 1, f_gstar2 = 0,
                Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)

# examples of estimates for ATE under f_gstar1 - f_gstar0:
tmleCom_wmT.bA.bY_SL.1a$ATE$estimates
tmleCom_wmT.bA.bY_SL.1a$ATE$vars

### 1.2 Community-level analysis with a pooled individual-level regression on outcome
## 1.2.1 glm using correctly specified Qform, misspecified hform.g0 and hform.gstar
tmleCom_Options(Qestimator = "glm__glm", gestimator = "glm__glm", maxNperBin = nrow(data))
tmleCom_wmT.bA.bY_Qcgm.1a <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), obs.wts = "equal.within.pop",
                community.step = "community_level", community.wts = "size.community", 
                communityID = "id", pooled.Q = FALSE, f_gstar1 = 1, f_gstar2 = 0,
                Qform = Qform.corr, hform.g0 = gform.mis, hform.gstar = gform.mis)
tmleCom_wmT.bA.bY_Qcgm.1a$ATE$estimates

#***************************************************************************************
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
#***************************************************************************************

#***************************************************************************************
# EXAMPLES OF Estimators:
#***************************************************************************************
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = nrow(dat_iidcontABinY))
tmleCom_res <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                             communityInd = NULL, community.step = NULL, f_gstar1 = f.gstar, Qform = Qform.corr, 
                             hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_res_same <- tmleCommunity(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                             f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_res_Alt <- tmleSingleStep(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                                  f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)

#***************************************************************************************
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = nrow(dat_iidcontABinY),
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), nbins = 10)
require(h2oEnsemble)
tmleCom_res <- tmleSingleStep(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                              f_gstar1 = f.gstar, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)    
