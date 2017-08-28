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
tmleCom_res <- tmleSingleStep(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                              f_gstar1 = f.gstar, Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)

#***************************************************************************************
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = nrow(dat_iidcontABinY),
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), nbins = 10)
require(h2oEnsemble)
tmleCom_res <- tmleSingleStep(data = dat_iidcontABinY, Ynode = "Y", Anodes = "A", Wnodes = c("W1", "W2", "W3", "W4"), 
                              f_gstar1 = f.gstar, Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)    
