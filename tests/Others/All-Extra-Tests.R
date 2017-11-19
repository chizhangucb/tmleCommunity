###############################################################################################
###############################################################################################
# tests_00_detailsExamples_tmleCommunity.R
###############################################################################################
###############################################################################################
# ---------------------------------------------------------------------------------
# TEST SET 0. SHOWING SOME DETAILS IN tmleCommunity 
library(data.table)
# ---------------------------------------------------------------------------------
# A is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
N <- nrow(indSample.iid.cA.cY)
head(indSample.iid.cA.cY)
psi0.Y <- mean(indSample.iid.cA.cY$Y)  # 0.291398
psi0.Ygstar <- mean(indSample.iid.cA.cY$Y.gstar)  # 0.316274
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
f.gstar <- define_f.gstar(shift = indSample.iid.cA.cY_list$shift.val, truncBD = indSample.iid.cA.cY_list$truncBD)

# ---------------------------------------------------------
# Step 1. Create an R6 object that stores and manages the input data, later passed on to estimation algorithm(s)
# ---------------------------------------------------------
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N, nbins = 10)
OData <- DatKeepClass$new(Odata = subset(indSample.iid.cA.cY, select=-Y), 
                          nodes = nodes[c("Anodes", "WEnodes")], norm.c.sVars = FALSE)
OData$nodes <- nodes
head(OData$dat.sVar)  # a subset data.frame of Odata that includes all variables in nodes
nobs <- OData$nobs
obsYvals <- indSample.iid.cA.cY[, nodes$Ynode]
ab <- range(obsYvals, na.rm=TRUE)
obsYvals <- (obsYvals-ab[1]) / diff(ab)
OData$addYnode(YnodeVals = obsYvals)
sum(OData$noNA.Ynodevals == obsYvals) == nobs  # TRUE
sum(OData$det.Y) == 0  # TRUE (No deterministic Y in input data)
OData$names.sVar  # "Y"  "A"  "W1" "W2" "W3" "W4" (names of all variables in input data )
OData$names.c.sVar  # "A"  "W3" "W4" (names of all continuous variables in input data)

# ---------------------------------------------------------
# Step 2. Defining and estimating outcome mechanism E(Y|A, E, W)
# ---------------------------------------------------------
Q.sVars <- tmleCommunity:::define_regform(NULL, Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "WEnodes")])
Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, predvars = Q.sVars$predvars, subset_vars = (!rep_len(FALSE, nobs)))
Qreg$estimator  # "speedglm__glm"
Qreg$pool_cont  # FALSE (Don't pool across all bins of a continuous exposure)
Qreg$S3class    # "RegressionClass" "R6"    

m.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData)
m.Q.init$getfit$coef
m.Q.init$predict(newdata = OData)
mean(m.Q.init$getprobA1 * diff(ab) + ab[1])  # 3.324724

# ---------------------------------------------------------
# Step 3. Defining and estimating treatment mechanism P(A|E, W) under g0
# ---------------------------------------------------------
h.g0.sVars <- tmleCommunity:::define_regform(NULL, Anodes.lst = nodes$Anodes, Wnodes.lst = nodes["WEnodes"])
sA_nms_g0 <- h.g0.sVars$outvars
subsets_expr <- lapply(sA_nms_g0, function(var) {var}) 
regclass.g0 <- RegressionClass$new(outvar = h.g0.sVars$outvars,
                                   predvars = h.g0.sVars$predvars,
                                   subset_vars = subsets_expr,
                                   outvar.class = OData$type.sVar[sA_nms_g0])
OData.g0 <- OData
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
genericmodels.g0$is.fitted  # FALSE
# genericmodels.g0$fit(data = OData.g0)
# h_gN <- genericmodels.g0$predictAeqa(newdata = OData); mean(h_gN)

# ---------------------------------------------- 
# Step 3.1 Details regarding A1 (the first exposure)
# ---------------------------------------------- 
genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()$`P(A|W).1` # genericmodels.g0$getPsAsW.models()[[1]]
genericmodels.g0.A1$intrvls
# Create a matrix of dummy bin indicators for continuous sVar 
OData.g0$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls, 
                       nbins = genericmodels.g0.A1$reg$nbins, bin.nms = genericmodels.g0.A1$bin_nms)  
bin.ind.mat <- OData.g0$dat.bin.sVar
colSums(bin.ind.mat, na.rm = T)
#  A_B.1  A_B.2  A_B.3  A_B.4  A_B.5  A_B.6  A_B.7  A_B.8  A_B.9 A_B.10 A_B.11 A_B.12 
#    0   1000   1000   1000   1000   1000   1000   1000   1000   1000   1000      0 
genericmodels.g0$fit(data = OData.g0)
genericmodels.g0.A1.B1 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[1]]
genericmodels.g0.A1.B2 <- genericmodels.g0$getPsAsW.models()[[1]]$getPsAsW.models()[[2]]
genericmodels.g0.A1.B2$getfit$coef  # the same as genericmodels.g0.A1.alt.B2$getfit$coef
#  Intercept         W1         W2         W3         W4 
# -0.6758011 -1.6536866 -0.7941279  0.5959984 -1.6902898

# Set regclass.g0 to a single univariate k_i regression for outcome regclass.g0$outvar[[k_i]]
k_i <- 1
reg_i <- regclass.g0$clone()
reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
genericmodels.g0.A1.alt <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData.g0)

# Define subset evaluation for new bins - used to initialize new GenericModel
subset_eval.A1 <- tmleCommunity:::def_regs_subset(genericmodels.g0.A1.alt)
subset_eval.A1$outvar.class  # a list of "binary", named as A_B.1, A_B.2, ..., A_B.12
subset_eval.A1$predvars  # "W1" "W2" "W3" "W4"
subset_eval.A1$subset_vars  #  a list of ("A_B.1", "A"), ("A_B.2", "A"), ..., ("A_B.12", "A"), named as A_B.1, A_B.2, ..., A_B.12

### 1) Fit ContinModel while saving dataset in each bin 
genericmodels.g0.A1.alt$fit(data = OData.g0, savespace = FALSE)
genericmodels.g0.A1.alt.B1 <- genericmodels.g0.A1.alt$getPsAsW.models()[[1]]
genericmodels.g0.A1.alt.B2 <- genericmodels.g0.A1.alt$getPsAsW.models()[[2]]
genericmodels.g0.A1.alt.B3 <- genericmodels.g0.A1.alt$getPsAsW.models()[[3]]
dim(genericmodels.g0.A1.alt.B2$getXmat)  # 10000 5
dim(genericmodels.g0.A1.alt.B3$getXmat)  # 9000 5
genericmodels.g0.A1.alt.B1$getfit$coef
#  Intercept            W1            W2            W3            W4 
# -2.656607e+01 -8.204045e-13 -1.361453e-13 -2.641735e-14 -2.408317e-12 
genericmodels.g0.A1.alt.B2$getfit$coef
#  Intercept         W1         W2         W3         W4 
# -0.6758011 -1.6536866 -0.7941279  0.5959984 -1.6902898 
length(genericmodels.g0.A1.alt.B3$getY)  # 9000
sum(genericmodels.g0.A1.alt.B3$subset_idx)  # 9000
head(genericmodels.g0.A1.alt.B3$getXmat)

### 2) Fit GenericModel directly from genericmodels.g0.A1.alt
genericmodels.g0.A1.alt.Generic <- GenericModel$new(reg = subset_eval.A1)
# The 1st bin with all Y = 0 (No observation falls into this protective bin)
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[1]]$fit(data = OData.g0, savespace = F)
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[1]]$getfit$coef  # the same as genericmodels.g0.A1.alt$getfit$coef
#  Intercept            W1            W2            W3            W4 
# -2.656607e+01 -8.204045e-13 -1.361453e-13 -2.641735e-14 -2.408317e-12 

# The 2nd bin (Real "first" bin)
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[2]]$fit(data = OData.g0, savespace = F)
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[2]]$getfit$coef  # The same result as genericmodels.g0.A1.alt$getfit$coef
#  Intercept         W1         W2         W3         W4 
# -0.6758011 -1.6536866 -0.7941279  0.5959984 -1.6902898 

# The 12th bin
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[12]]$fit(data = OData.g0, savespace = F)
genericmodels.g0.A1.alt.Generic$getPsAsW.models()[[12]]$getfit$coef  # NA NA NA NA NA

# ---------------------------------------------- 
# Step 3.2 Details regarding get.sVar.bw
# ---------------------------------------------- 
# Reinitiate genericmodels.g0 to avoild stopping when is.fitted = TRUE & overwrite = FALSE
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)  
genericmodels.g0$fit(data = OData.g0)
genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()[[1]]
OData.new <- OData
OData.new$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls, 
                        nbins = genericmodels.g0.A1$nbins, bin.nms = genericmodels.g0.A1$bin_nms)
# Create a vector of ordinal (categorical) vars out of cont. sVar vector
ord.sVar <- tmleCommunity:::make.ordinal(x = OData.new$get.sVar(genericmodels.g0.A1$outvar), 
                                         intervals = genericmodels.g0.A1$intrvls)
head(ord.sVar, 10)
table(ord.sVar)
#    2    3    4    5    6    7    8    9   10   11 
# 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
sum(ord.sVar == 1)   # 0  Indicate NO data point falls beyond the minimum bound
sum(ord.sVar == 12)  # 0  Indicate NO data point falls beyond the maximum bound
head(OData.new$get.sVar(genericmodels.g0.A1$outvar), 10)  # Get outcome vector
intrvls.width <- diff(genericmodels.g0.A1$intrvls)
bws <- OData.new$get.sVar.bw(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls)
genericmodels.g0.A1$bin_weights <- bin_weights <- (1 / bws)

# ---------------------------------------------- 
# Step 3.3 Details regarding predictAeqa
# ---------------------------------------------- 
# Reclone regclass.g0 to avoid unexpected problem
k_i <- 1
reg_i <- regclass.g0$clone()
reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
datX_mat <- OData.new$get.dat.sVar(rowsubset = TRUE, covars = (genericmodels.g0.A1$bin_nms))
pooled_bin_name <- OData.new$pooled.bin.nm.sVar(reg_i$outvar)  # "A_allB.j"
genericmodels.g0.A1.alt <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData.new)
# Convert Bin matrix for continuous exposure into long format data.table
BinsDat_long <- tmleCommunity:::binirized.to.DTlong(
  BinsDat_wide = datX_mat, binID_seq = (1L:reg_i$nbins), ID = as.integer(1:nrow(datX_mat)),
  bin_names = genericmodels.g0.A1$bin_names, pooled_bin_name = pooled_bin_name, name.sVar = reg_i$outvar)
sVar_melt_DT <- tmleCommunity:::join.Xmat(X_mat = OData.new$get.dat.sVar(TRUE, covars = reg_i$predvars), 
                                          sVar_melt_DT = BinsDat_long, ID = as.integer(1:nrow(datX_mat)))
X_mat <- sVar_melt_DT[,c("bin_ID", reg_i$predvars), with=FALSE][, c("Intercept") := 1] # select bin_ID + predictors, add intercept column
setcolorder(X_mat, c("Intercept", "bin_ID", reg_i$predvars)) # re-order columns by reference (no copy)
dim(X_mat)  # 65000  6
ID <- sVar_melt_DT[["ID"]]
Y_vals <- sVar_melt_DT[, pooled_bin_name, with = FALSE][[1]] 

# ---------------------------------------------- 
# Step 3.4 Details regarding predictAeqa
# ---------------------------------------------- 
genericmodels.g0.A1.B2 <- genericmodels.g0.A1$getPsAsW.models()[[2]]
genericmodels.g0.A1.B2$newdata(newdata = OData, getoutvar = TRUE)
genericmodels.g0.A1.B2$getfit$fitfunname  # "speedglm"
probA1 <- tmleCommunity:::predict_single_reg(self = genericmodels.g0.A1.B2)
indA <- OData.new$get.outvar(genericmodels.g0.A1.B2$subset_idx, var = genericmodels.g0.A1.B2$outvar)
probAeqa <- rep.int(1L, nobs) 
probA1 <- probA1[genericmodels.g0.A1.B2$subset_idx]
probAeqa[genericmodels.g0.A1.B2$subset_idx] <- probA1^(indA) * (1 - probA1)^(1L - indA)
mean(probAeqa)  # 0.7192802

## Test if the cumulative probability * bin weights will be the same result as $getcumprodAeqa() provides
genericmodels.g0.A1$predictAeqa(newdata = OData, wipeProb = FALSE)
mean(genericmodels.g0.A1$getcumprodAeqa())  # 0.2345022
cumprodAeqa <- 1
for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
  cumprodAeqa <- cumprodAeqa * genericmodels.g0.A1$getPsAsW.models()[[i]]$getprobAeqa
}
mean(cumprodAeqa * bin_weights) == mean(genericmodels.g0.A1$getcumprodAeqa())  # TRUE

genericmodels.g0.A1.B2 <- genericmodels.g0.A1$getPsAsW.models()[[2]]
head(genericmodels.g0.A1.B2$getprobA1); head(genericmodels.g0.A1.B2$getprobAeqa)
mean(genericmodels.g0.A1.B2$getprobAeqa)  # 0.7192802
sapply(1:length(genericmodels.g0.A1$getPsAsW.models()), FUN = function(x) {
  mean(genericmodels.g0.A1$getPsAsW.models()[[x]]$getprobAeqa)
})  # Provide every cumulative probability for each bin
#  1.0000000 0.7192802 0.7213857 0.7501958 0.8099905 1.0000000 1.0000000


# ---------------------------------------------------------
# Step 4. Defining and estimating treatment mechanism P(A|E, W) under gstar
# ---------------------------------------------------------
h.gstar.sVars <- h.g0.sVars
sA_nms_gstar <- h.gstar.sVars$outvars
regclass.gstar <- RegressionClass$new(outvar = sA_nms_gstar,
                                      predvars = h.gstar.sVars$predvars,
                                      subset_vars = subsets_expr,
                                      outvar.class = OData$type.sVar[sA_nms_gstar])
OData.gstar <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes, norm.c.sVars = FALSE)
OData.gstar$make.dat.sVar(p = 1, f.g_fun = f.gstar)

# -------------------------------------------------------------------------------------------
# estimating h_g0 and h_gstar
# -------------------------------------------------------------------------------------------
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
genericmodels.g0$fit(data = OData.g0)
h_gN <- genericmodels.g0$predictAeqa(newdata = OData)

genericmodels.gstar <- GenericModel$new(reg = regclass.gstar, DatKeepClass.g0 = OData.g0)
genericmodels.gstar$fit(data = OData.gstar)
h_gstar <- genericmodels.gstar$predictAeqa(newdata = OData)


###############################################################################################
###############################################################################################
# bootstrap_tests.R
###############################################################################################
###############################################################################################
######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
truncBD <- 5
shift.val <- 1

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    # print("shift.const: " %+% shift.const)
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
                           f.gstar, lbound = 0.005, working.model = TRUE, Qform = NULL, gform = NULL, rndseed = NULL, verbose = F) {
  
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
    data <- try(get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = NULL, truncBD = truncBD, shift.val = shift.val, verbose = verbose,
                                  is.Y.bin = TRUE, working.model = working.model, timevarying = TRUE, n.ind.fix = FALSE))
    if (inherits(data, "try-error")) { error_dat <- error_dat + 1; next } # skip iterations with errors
    
    tmleCom_Options(Qestimator = Qestimator, gestimator = gestimator, maxNperBin = nrow(data), nbins = nbins)
    
    # TMLE-Ia
    tmleCom_1a_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "community_level", working.model = FALSE, f_gstar1 = f.gstar,
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound, verbose = verbose))
    if (inherits(tmleCom_1a_size, "try-error")) { error_1a <- error_1a + 1; next } # skip iterations with errors
    message("################### Finishing the TMLE-Ia ###################\n")
    
    # TMLE-Ib
    tmleCom_1b_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "individual_level", working.model = FALSE,  f_gstar1 = f.gstar, 
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound, verbose = verbose))
    if (inherits(tmleCom_1b_size, "try-error")) { error_1b <- error_1b + 1; next } # skip iterations with errors
    message("################### Finishing the TMLE-Ib ###################\n")
    
    # TMLE-II
    tmleCom_ii_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                         obs.wts = "equal.within.community", community.wts = "size.community", communityID = "id",
                                         community.step = "individual_level", working.model = TRUE, f_gstar1 = f.gstar,
                                         Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound, verbose = verbose))
    if (inherits(tmleCom_ii_size, "try-error")) { error_2 <- error_2 + 1; next } # skip iterations with errors
    message("################### Finishing the TMLE-II ###################\n")
    
    # TMLE-Per
    tmleCom_per_size <- try(tmleCommunity(data = data, Ynode = "Y", Anodes = "A", WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                          obs.wts = "equal.within.community", community.wts = "equal.community", communityID = "id",
                                          community.step = "perCommunity", working.model = TRUE, f_gstar1 = f.gstar,
                                          Qform = Qform, hform.g0 = gform, hform.gstar = gform, lbound = lbound, verbose = verbose))
    if (inherits(tmleCom_per_size, "try-error")) { error_Per <- error_Per + 1; next } # skip iterations with errors
    message("################### Finishing the TMLE-Stratifed ###################\n")
    
    est_tmle1a <- rbind(est_tmle1a, tmleCom_1a_size$EY_gstar1$estimates[, 1])
    var_tmle1a <- rbind(var_tmle1a, tmleCom_1a_size$EY_gstar1$vars[, 1])
    CI_tmle1a <- rbind(CI_tmle1a, c(tmleCom_1a_size$EY_gstar1$CIs[, 1], tmleCom_1a_size$EY_gstar1$CIs[, 2]))
    pval_tmle1a <- rbind(pval_tmle1a, tmleCom_1a_size$EY_gstar1$pval[, 1])
    
    est_tmle1b <- rbind(est_tmle1b, tmleCom_1b_size$EY_gstar1$estimates[, 1])
    var_tmle1b <- rbind(var_tmle1b, tmleCom_1b_size$EY_gstar1$vars[, 1])
    CI_tmle1b <- rbind(CI_tmle1b, c(tmleCom_1b_size$EY_gstar1$CIs[, 1], tmleCom_1b_size$EY_gstar1$CIs[, 2]))
    pval_tmle1b <- rbind(pval_tmle1b, tmleCom_1b_size$EY_gstar1$pval[, 1])
    
    est_tmle2 <- rbind(est_tmle2, tmleCom_ii_size$EY_gstar1$estimates[, 1])
    var_tmle2 <- rbind(var_tmle2, tmleCom_ii_size$EY_gstar1$vars[, 1])
    CI_tmle2 <- rbind(CI_tmle2, c(tmleCom_ii_size$EY_gstar1$CIs[, 1], tmleCom_ii_size$EY_gstar1$CIs[, 2]))
    pval_tmle2 <- rbind(pval_tmle2, tmleCom_ii_size$EY_gstar1$pval[, 1])
    
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


#######################################################
############# Whne working model fails ################
#######################################################
## Test 1.1.1 Correctly specified Qform & gform
## Test 1.1.1.1 The number of communities = 30 & the number of individual per community = 30
tmle.wmF.Qc.gc.J30.n30.T <- 
  system.time(tmle.wmF.Qc.gc.J30.n30 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 30, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                       Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                       lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J30.n30, tmle.wmF.Qc.gc.J30.n30.T, file = "tmle.wmF.Qc.gc.J30.n30.Rda")

## Test 1.1.1.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmF.Qc.gc.J300.n300.T <- 
  system.time(tmle.wmF.Qc.gc.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 200, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J300.n300, tmle.wmF.Qc.gc.J300.n300.T, file = "tmle.wmF.Qc.gc.J300.n300.Rda")

## Test 1.1.1.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmF.Qc.gc.J500.n300.T <-
  system.time(tmle.wmF.Qc.gc.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J500.n300, tmle.wmF.Qc.gc.J500.n300.T, file = "tmle.wmF.Qc.gc.J500.n300.Rda")

## Test 1.1.1.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmF.Qc.gc.J1000.n300.T <- 
  system.time(tmle.wmF.Qc.gc.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J1000.n300, tmle.wmF.Qc.gc.J1000.n300.T, file = "tmle.wmF.Qc.gc.J1000.n300.Rda")

## Test 1.1.1.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmF.Qc.gc.J300.n500.T <- 
  system.time(tmle.wmF.Qc.gc.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J300.n500, tmle.wmF.Qc.gc.J300.n500.T, file = "tmle.wmF.Qc.gc.J300.n500.Rda")

## Test 1.1.1.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmF.Qc.gc.J500.n500.T <-
  system.time(tmle.wmF.Qc.gc.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qc.gc.J500.n500, tmle.wmF.Qc.gc.J500.n500.T, file = "tmle.wmF.Qc.gc.J500.n500.Rda")

## Test 1.1.2 Correctly specified Qform & Misspecified gform
## Test 1.1.2.1 The number of communities = 30 & the number of individual per community = 300
tmle.wmF.Qc.gm.J30.n300.T <- 
  system.time(tmle.wmF.Qc.gm.J30.n300 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                        Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                        lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J30.n300, tmle.wmF.Qc.gm.J30.n300.T, file = "tmle.wmF.Qc.gm.J30.n300.Rda")

## Test 1.1.2.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmF.Qc.gm.J300.n300.T <-
  system.time(tmle.wmF.Qc.gm.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J300.n300,tmle.wmF.Qc.gm.J300.n300.T, file = "tmle.wmF.Qc.gm.J300.n300.Rda")

## Test 1.1.2.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmF.Qc.gm.J500.n300.T <- 
  system.time(tmle.wmF.Qc.gm.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J500.n300, tmle.wmF.Qc.gm.J500.n300.T, file = "tmle.wmF.Qc.gm.J500.n300.Rda")

## Test 1.1.2.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmF.Qc.gm.J1000.n300.T <- 
  system.time(tmle.wmF.Qc.gm.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J1000.n300, tmle.wmF.Qc.gm.J1000.n300.T, file = "tmle.wmF.Qc.gm.J1000.n300.Rda")

## Test 1.1.2.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmF.Qc.gm.J300.n500.T <- 
  system.time(tmle.wmF.Qc.gm.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J300.n500, tmle.wmF.Qc.gm.J300.n500.T, file = "tmle.wmF.Qc.gm.J300.n500.Rda")

## Test 1.1.2.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmF.Qc.gm.J500.n500.T <-
  system.time(tmle.wmF.Qc.gm.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.corr, gform = gform.mis, rndseed = 1))
save(tmle.wmF.Qc.gm.J500.n500, tmle.wmF.Qc.gm.J500.n500.T, file = "tmle.wmF.Qc.gm.J500.n500.Rda")

## Test 1.1.3 Misspecified Qform & Correctly specified gform
## Test 1.1.3.1 The number of communities = 30 & the number of individual per community = 300
tmle.wmF.Qm.gc.J30.n300.T <- 
  system.time(tmle.wmF.Qm.gc.J30.n300 <- bootstrap.TMLE(nBoot = 200, J = 30, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                        Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                        lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J30.n300, tmle.wmF.Qm.gc.J30.n300.T, file = "tmle.wmF.Qm.gc.J30.n300.Rda")

## Test 1.1.3.2 The number of communities = 300 & the number of individual per community = 300
tmle.wmF.Qm.gc.J300.n300.T <- 
  system.time(tmle.wmF.Qm.gc.J300.n300 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J300.n300, tmle.wmF.Qm.gc.J300.n300.T, file = "tmle.wmF.Qm.gc.J300.n300.Rda")

## Test 1.1.3.3 The number of communities = 500 & the number of individual per community = 300
tmle.wmF.Qm.gc.J500.n300.T <- 
  system.time(tmle.wmF.Qm.gc.J500.n300 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J500.n300, tmle.wmF.Qm.gc.J500.n300.T, file = "tmle.wmF.Qm.gc.J500.n300.Rda")

## Test 1.1.3.4 The number of communities = 1000 & the number of individual per community = 300
tmle.wmF.Qm.gc.J1000.n300.T <- 
  system.time(tmle.wmF.Qm.gc.J1000.n300 <- bootstrap.TMLE(nBoot = 200, J = 1000, n.ind = 300, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                          Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                          lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J1000.n300, tmle.wmF.Qm.gc.J1000.n300.T, file = "tmle.wmF.Qm.gc.J1000.n300.Rda")

## Test 1.1.3.5 The number of communities = 300 & the number of individual per community = 500
tmle.wmF.Qm.gc.J300.n500.T <- 
  system.time(tmle.wmF.Qm.gc.J300.n500 <- bootstrap.TMLE(nBoot = 200, J = 300, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J300.n500, tmle.wmF.Qm.gc.J300.n500.T, file = "tmle.wmF.Qm.gc.J300.n500.Rda")

## Test 1.1.3.6 The number of communities = 500 & the number of individual per community = 500
tmle.wmF.Qm.gc.J500.n500.T <- 
  system.time(tmle.wmF.Qm.gc.J500.n500 <- bootstrap.TMLE(nBoot = 200, J = 500, n.ind = 500, truncBD = truncBD, shift.val = shift.val, nbins = 5,
                                                         Qestimator = "speedglm__glm", gestimator = "speedglm__glm", f.gstar = f.gstar, 
                                                         lbound = 0.025, working.model = FALSE, Qform = Qform.mis, gform = gform.corr, rndseed = 1))
save(tmle.wmF.Qm.gc.J500.n500, tmle.wmF.Qm.gc.J500.n500.T, file = "tmle.wmF.Qm.gc.J500.n500.Rda")


###############################################################################################
###############################################################################################
# tmleComTests.R
###############################################################################################
###############################################################################################
load("Case I/EDA/data.Rda")
YsName <- c("premature_est", "zscore_bwgtgest")
AsName <- c("viol_r_tri1", "viol_r_tri2")
covEfix.Name <- c("pct_foreign_born", "pct_poverty", "pct_black", "pct_hisp", "pct_renters")
covEextra.Name <- c("umemp_a", "maxt_a", "pcp_a")
covEName <- c(covEfix.Name, covEextra.Name)
covWName <- c(c(paste0("mrace", 1:5), paste0("meduc_cat", 1:4), paste0("month", 1:11), 
                paste0("yr", 2007:2011)), "age_mom", "parity", "pay_gov_ins", "foreign_born")

create_f.gstar <- function(shiftPercent, truncBD = NULL, truncFUN = NULL, covNames) {
  shift.percent <- shiftPercent
  if (is.function(truncFUN)) { 
    truncFUN <- truncFUN
  } else if (all(is.number(truncBD))) {
    trunc.bound <- truncBD
  } else {
    stop("Either truncBD or truncFUN should be defined," %+%
           "where truncBD is a list of numbers and truncFUN is a function")
  }
  trtsName <- covNames
  f.gstar <- function(data, ...) {
    if (is.null(truncBD)) trunc.bound <- c(truncFUN(data[, trtsName[1]]), truncFUN(data[, trtsName[2]]))
    print("shift percentage: " %+% shift.percent)
    print("truncated bound for " %+% trtsName[1] %+% ": " %+% trunc.bound[1])
    print("truncated bound for " %+% trtsName[2] %+% ": " %+% trunc.bound[1])
    untrunc.A1 <- data[, trtsName[1]] * shift.percent
    untrunc.A2 <- data[, trtsName[2]] * shift.percent
    trunc.A1 <- base:::ifelse(trunc.bound[1] > untrunc.A1, trunc.bound[1], untrunc.A1)
    trunc.A2 <- base:::ifelse(trunc.bound[2] > untrunc.A2, trunc.bound[2], untrunc.A2)
    return(cbind(trunc.A1, trunc.A2))
  }
  return(f.gstar)
}
f.gstar <- create_f.gstar(shiftPercent = 0.8, truncFUN = min, covNames = AsName)

data <- subdata
gvars$verbose <- TRUE

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
tmleCom_Options(maxNperBin = nrow(data))
tmleCommunity_res.glm1 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = f.gstar, Qform = NULL, n_MCsims = 1)

tmleCommunity_res.glm1$EY_gstar1$estimates
#        estimate
# tmle   0.07298817
# h_iptw 0.07444963
# gcomp  0.07021494

tmleCommunity_res.glm2 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCommunity_res.glm2$EY_gstar1$estimates
#        estimate
# tmle    0.07016
# h_iptw  0.07016
# gcomp   0.07016
tmleCommunity_res.glm2$EY_gstar2$estimates
#        estimate
# tmle   0.07298817
# h_iptw 0.07444963
# gcomp  0.07021494
tmleCommunity_res.glm2$ATE$estimates
#             estimate
# tmle   -2.828174e-03
# h_iptw -4.289633e-03
# gcomp  -5.494195e-05


######################################### 
## Test 1.2 h2o
#########################################
tmleCom_Options(estimator = "h2o__ensemble", h2olearner =  c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), maxNperBin = nrow(data))
tmleCom_res.h2o.1 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCom_res.h2o.1$EY_gstar1$estimates
#          estimate
# tmle   0.07016000
# h_iptw 0.07016000
# gcomp  0.07070904
tmleCom_res.h2o.1$EY_gstar2$estimates
#           estimate
# tmle   0.071036598
# h_iptw 0.009401554
# gcomp  0.070604479
tmleCom_res.h2o.1$ATE$estimates
#             estimate
# tmle   -0.0008765975
# h_iptw  0.0607584455
# gcomp   0.0001045629

tmleCom_Options(estimator = "h2o__ensemble", maxNperBin = nrow(data),
                h2olearner =  as.vector(paste('h2o.glm.', c(0, 0.50, 1), sep = "")))
tmleCom_res.h2o.2 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                   Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCom_res.h2o.2$EY_gstar1$estimates
#          estimate
# tmle   0.0701600
# h_iptw 0.07016000
# gcomp  0.07011594
tmleCom_res.h2o.2$EY_gstar2$estimates
#          estimate
# tmle   0.07096359
# h_iptw 0.08344168
# gcomp  0.07021925
tmleCom_res.h2o.2$ATE$estimates
#             estimate
# tmle   -0.0008035862
# h_iptw -0.0132816815
# gcomp  -0.0001033166
