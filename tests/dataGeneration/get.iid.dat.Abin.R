# ---------------------------------------------------------------------------------
# DATA SET 2. TESTS FOR FITTING BINARY POINT EXPOSURE A IN IID DATA
# ---------------------------------------------------------------------------------
# Fitting binary exposure by logistic regression
# individual exposure is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
# ---------------------------------------------------------------------------------

get.iid.dat.Abin <- function(ndata = 100000, rndseed = NULL, is.Y.bin = TRUE) {
  require(simcausal)
  D <- DAG.empty()
  D <- D + 
    node("W1", distr = "rbern", prob = 0.5) +
    node("W2", distr = "rbern", prob = 0.3) +
    node("W3", distr = "rnorm", mean = 0, sd = 0.3) +
    node("W4", distr = "runif", min = 0, max = 1) +
    node("W3W4", distr = "rconst", const = W3 * W4) +
    node("A", distr = "rbern", prob = plogis(0.86 * W1 + W2 + 1.32 * W3 - W4 - 0.45 *W3W4 + 0.1)) +
    node("W2A", distr = "rconst", const = W2 * A)
  
  if (is.Y.bin) {  # Create binary Y
    D <- D +  
      node("Y", distr = "rbern", prob = plogis(2.8 * A + 2 * W1 + 1.5 * W2 + 0.8 * W2A  + 0.5 * W3 - 3 * W4 - 3.5)) + 
      node("Y1", distr = "rbern", prob = plogis(2.8 * 1 + 2 * W1 + 1.5 * W2 + 0.8 * W2 * 1  + 0.5 * W3 - 3 * W4 - 3.5)) +
      node("Y0", distr = "rbern", prob = plogis(2.8 * 0 + 2 * W1 + 1.5 * W2 + 0.8 * W2 * 0  + 0.5 * W3 - 3 * W4 - 3.5))
  } else {  # Create continuous Y
    D <- D +  
      node("Y", distr = "rnorm", mean = (2.8 * A + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) + 
      node("Y1", distr = "rnorm", mean = (2.8 * 1 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) +
      node("Y0", distr = "rnorm", mean = (2.8 * 0 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1)
  }
  
  D <- set.DAG(D)
  Odata <- sim(D, n = ndata, rndseed = rndseed)
  psi0.Y <- mean(Odata$Y1) - mean(Odata$Y0)
  print("psi0.Y: " %+% psi0.Y)  
  return(list(psi0.Y = psi0.Y, Odata = Odata))
}

ndata <- 10000
rndseed <- 12345

#### Data 1. One binary, individual-level A with binary Y
indPop.ind.bA.bY <- get.iid.dat.Abin(ndata = 1000000, rndseed = rndseed, is.Y.bin = TRUE)$Odata
psi0.Y <- mean(indPop.ind.bA.bY$Y1) - mean(indPop.ind.bA.bY$Y0) # 0.348242
indSample.ind.bA.bY <- get.iid.dat.Abin(ndata = ndata, rndseed = rndseed, is.Y.bin = TRUE)$Odata
indSample.ind.bA.bY <- indSample.ind.bA.bY[, c("W1", "W2", "W3", "W4", "A", "Y")]
indSample.ind.bA.bY_list <- list(indSample.ind.bA.bY = indSample.ind.bA.bY, rndseed = rndseed, psi0.Y = psi0.Y,
                                 psi0.Y1 = mean(indPop.ind.bA.bY$Y1), psi0.Y0 = mean(indPop.ind.bA.bY$Y0))
save(indSample.ind.bA.bY_list, file="indSample.ind.bA.bY_list.Rda")

#### Data 2. One binary, individual-level A with continuous Y
indPop.ind.bA.cY <- get.iid.dat.Abin(ndata = 1000000, rndseed = rndseed, is.Y.bin = FALSE)$Odata
psi0.Y <- mean(indPop.ind.bA.cY$Y1) - mean(indPop.ind.bA.cY$Y0) # 1.800267
indSample.ind.bA.cY <- get.iid.dat.Abin(ndata = ndata, rndseed = rndseed, is.Y.bin = FALSE)$Odata
indSample.ind.bA.cY <- indSample.ind.bA.cY[, c("W1", "W2", "W3", "W4", "A", "Y")]
indSample.ind.bA.cY_list <- list(indSample.ind.bA.cY = indSample.ind.bA.cY, rndseed = rndseed, psi0.Y = psi0.Y,
                                 psi0.Y1 = mean(indPop.ind.bA.cY$Y1), psi0.Y0 = mean(indPop.ind.bA.cY$Y0))
save(indSample.ind.bA.cY_list, file="indSample.ind.bA.cY_list.Rda")


################################################
#### FUN 2. Rare outcome with binary Y (Case-control study)
################################################
get.iid.dat.Abin2 <- function(ndata = 100000, rndseed = NULL) {
  require(simcausal)
  D <- DAG.empty()
  D <- D + 
    node("W1", distr = "runif", min = 0, max = 1) +
    node("W2", distr = "rnorm", mean = 0, sd = 0.3) +
    node("W3", distr = "rbern", prob = 0.5) +
    node("W4", distr = "rbern", prob = 0.5) +
    node("A", distr = "rbern", prob = plogis(W1 + W2 + 2 * W3 + W4 - 0.8)) +
    node("W2A", distr = "rconst", const = W2 * A) +
    node("Y", distr = "rbern", prob = plogis(1.8 * A + 2 * W1 + W2 + 0.8 * W2A  + 0.5 * W3 - 4 * W4 - 3.5) / 12) + 
    node("Y1", distr = "rbern", prob = plogis(1.8 * 1 + 2 * W1 + W2 + 0.8 * W2 * 1  + 0.5 * W3 - 4 * W4 - 3.5) / 12) +
    node("Y0", distr = "rbern", prob = plogis(1.8 * 0 + 2 * W1 + W2 + 0.8 * W2 * 0  + 0.5 * W3 - 4 * W4 - 3.5) / 12)
  
  D <- set.DAG(D)
  Odata <- sim(D, n = ndata, rndseed = rndseed)
  print("mean(Odata$Y): " %+% mean(Odata$Y))
  psi0.Y <- mean(Odata$Y1) - mean(Odata$Y0)
  print("psi0.Y: " %+% psi0.Y)  
  
  return(list(psi0.Y = psi0.Y, Odata = Odata))
}

drawSamples <- function(nCa, nCo, Cohort, rndseed = NULL, replace = FALSE) {
  ## Description: Draw case-control samples 
  ## nCa: the number of cases 
  ## nCo: the number of controls 
  ## Cohort: a population data frame 
  
  q0 <- mean(Cohort$Y)
  n <- nCa + nCo
  J <- nCo / nCa  # ratio of control to case
  
  # Draw the case-control sample
  set.seed(rndseed)
  cases.ind <- sample(which(Cohort$Y == 1), size = nCa, replace = replace)
  controls.ind <- sample(which(Cohort$Y == 0), size = nCo, replace = replace)
  samples <- Cohort[c(cases.ind, controls.ind), ]
  
  # Assign weights to individuals 
  obs.wt <- rep(NA, n)
  obs.wt[samples$Y == 1] <- q0
  obs.wt[samples$Y == 0] <- (1 - q0) / J
  return(list(samples = samples, obs.wt = obs.wt, J = J))
}

rndseed <- 12345
cohort_iidBinABinY.rare <- get.iid.dat.Abin2(ndata = 1000000, rndseed = rndseed)$Odata
q0 <- mean(cohort_iidBinABinY.rare$Y)  # 0.013579 Rare outcome!!
psi0.Y <- mean(cohort_iidBinABinY.rare$Y1) - mean(cohort_iidBinABinY.rare$Y0) # 0.012662

# ------------------------------
# Sample set 1. J = 1
# ------------------------------
drawSamples.J1 <- drawSamples(nCa = 1000, nCo = 1000, rndseed = rndseed, replace = FALSE,
                              Cohort = cohort_iidBinABinY.rare[, c("W1", "W2", "W3", "W4", "A", "Y")])
dat_iidBinAY.rareJ1 <- drawSamples.J1$samples
obs.wt.J1 <- drawSamples.J1$obs.wt
rareSamples_iidBinAY.J1 <- list(dat_iidBinAY.rareJ1 = dat_iidBinAY.rareJ1, obs.wt.J1 = obs.wt.J1,
                                J = drawSamples.J1$J, q0 = q0, psi0.Y = psi0.Y, rndseed = rndseed, 
                                call = list(Cohort = "get.iid.dat.Abin2(ndata = 1000000, rndseed = 12345)$Odata",
                                            Sample = "drawSamples(nCa = 1000, nCo = 1000, rndseed = 12345, Cohort = Cohort)"))
save(rareSamples_iidBinAY.J1, file = "rareSamples_iidBinAY.J1.Rda")

# ------------------------------
# Sample set 1. J = 2
# ------------------------------
drawSamples.J2 <- drawSamples(nCa = 1000, nCo = 2000, rndseed = rndseed, replace = FALSE,
                              Cohort = cohort_iidBinABinY.rare[, c("W1", "W2", "W3", "W4", "A", "Y")])
dat_iidBinAY.rareJ2 <- drawSamples.J2$samples
obs.wt.J2 <- drawSamples.J2$obs.wt
rareSamples_iidBinAY.J2 <- list(dat_iidBinAY.rareJ2 = dat_iidBinAY.rareJ2, obs.wt.J2 = obs.wt.J2,
                                J = drawSamples.J2$J, q0 = q0, psi0.Y = psi0.Y, rndseed = rndseed, 
                                call = list(Cohort = "get.iid.dat.Abin2(ndata = 1000000, rndseed = 12345)$Odata",
                                            Sample = "drawSamples(nCa = 1000, nCo = 2000, rndseed = 12345, Cohort = Cohort)"))
save(rareSamples_iidBinAY.J2, file = "rareSamples_iidBinAY.J2.Rda")


################################################
#### FUN 3. One Binary A with Continuous Y
################################################
get.iid.dat.Abin3 <- function(ndata = 100000, rndseed = NULL) {
  require(simcausal)
  D <- DAG.empty()
    D <- D + 
    node("W1", distr = "rbern", prob = 0.5) +
    node("W2", distr = "rbern", prob = 0.3) +
    node("W3", distr = "rnorm", mean = 0, sd = 0.3) +
    node("W4", distr = "runif", min = 0, max = 1) +
    node("A", distr = "rbern", prob = plogis(0.92 * W1 + W2 + 1.32 * W3 - W4  + 0.1)) +
    node("Y", distr = "rnorm", mean = (2.8 * A + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) + 
    node("Y1", distr = "rnorm", mean = (2.8 * 1 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1) +
    node("Y0", distr = "rnorm", mean = (2.8 * 0 + 2 * W1 + 1.5 * W2 + 0.5 * W3 - 3 * W4 - 3.5), sd = 1)
  
  D <- set.DAG(D)
  Odata <- sim(D, n = ndata, rndseed = rndseed)
  psi0.Y <- mean(Odata$Y1) - mean(Odata$Y0)
  print("psi0.Y: " %+% psi0.Y)  
  
  return(list(psi0.Y = psi0.Y, Odata = Odata))
}

ndata <- 10000
rndseed <- 12345
popDat <- get.iid.dat.Abin3(ndata = 1000000, rndseed = rndseed, is.Y.bin = TRUE)$Odata
psi0.Y <- mean(popDat$Y1) - mean(popDat$Y0) # 0.348242

dat_iidBinAContY <- get.iid.dat.Abin3(ndata = ndata, rndseed = rndseed)$Odata
dat_iidBinAContY <- dat_iidBinAContY[, c("W1", "W2", "W3", "W4", "A", "Y")]
sampleDat_iidBinAContY <- list(dat_iidBinAContY = dat_iidBinAContY, psi0.Y = psi0.Y, rndseed = rndseed,
                               call = "get.iid.dat.Abin3(ndata = 10000, rndseed = 12345)$Odata")
save(sampleDat_iidBinAContY, file="sampleDat_iidBinAContY.Rda")


