# ---------------------------------------------------------------------------------
# DATA SET 1. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN IID DATA (No Cluster)
# ---------------------------------------------------------------------------------
# Fitting continuous exposure by binning, conditional on covariates (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (W1, W2, W3, W4), sd = 1;
# ---------------------------------------------------------------------------------

get.iid.dat.Acont <- function(ndata = 100000, rndseed = NULL, truncBD = 10, shift.val = 2, is.Y.bin = T) {
  require(simcausal)
  D <- DAG.empty()
  D <- D + 
    node("W1", distr = "rbern", prob = 0.5) +
    node("W2", distr = "rbern", prob = 0.3) +
    node("W3", distr = "rnorm", mean = 0, sd = 0.25) +
    node("W4", distr = "runif", min = 0, max = 1) +
    node("A.mu", distr = "rconst", const = (0.86 * W1 + 0.41 * W2 - 0.34 * W3 + 0.93 * W4)) +
    node("A", distr = "rnorm", mean = A.mu, sd = 1) +
    node("shift", distr = "rconst", const = .(shift.val)) +
    node("untrunc.A.gstar", distr = "rconst", const = A + shift) +
    node("trunc.c", distr = "rconst", const = .(truncBD))
    
  if (is.Y.bin) {  # Create binary Y
    D <- D + 
      node("r.new.A", distr = "rconst", const = exp(0.8 * shift * (untrunc.A.gstar - A.mu - shift / 3))) +
      node("trunc.A.gstar",  distr = "rconst", const = ifelse(r.new.A > trunc.c, A, untrunc.A.gstar)) + 
      node("Y", distr = "rbern", prob = plogis(-0.22 + 0.12 * A - 0.92 * W1 - 0.36 * W2 + 0.12 * W3 - 0.53 * W4)) +
      node("Y.gstar", distr = "rbern", prob = plogis(-0.22 + 0.12 *trunc.A.gstar - 0.92 * W1 - 0.36 * W2 + 0.12 * W3 - 0.53 * W4))
  } else {  # Create continuous Y
    D <- D + 
      node("r.new.A", distr = "rconst", const = exp(0.5 * shift * (untrunc.A.gstar - A.mu - shift / 2))) +
      node("trunc.A.gstar", distr = "rconst", const = ifelse(r.new.A > trunc.c, A, untrunc.A.gstar)) + 
      node("Y", distr = "rnorm", mean = (3.63 + 0.11 * A - 0.52 * W1 - 0.36 * W2 + 0.12 * W3 - 0.13 * W4), sd = 1) +
      node("Y.gstar", distr = "rnorm", mean = (3.63 + 0.11 * trunc.A.gstar - 0.52 * W1 - 0.36 * W2 + 0.12 * W3 - 0.13 * W4), sd = 1)
  }
 
  D <- set.DAG(D)
  Odata <- sim(D, n = ndata, rndseed = rndseed)
  # head(Odata, 50)
  psi0.Y <- mean(Odata$Y)
  print("psi0.Y: " %+% psi0.Y)  
  psi0.Ygstar <- mean(Odata$Y.gstar)
  print("psi0.Ygstar: " %+% psi0.Ygstar)
  
  return(list(psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar, Odata = Odata))
}

ndata <- 10000
rndseed <- 12345
truncBD <- 10

#### Data 1. One continuous A with Binary Y
shift.val <- 1
indPop.iid.cA.bY <- get.iid.dat.Acont(ndata = 1000000, rndseed = rndseed, truncBD = truncBD, shift.val = shift.val, is.Y.bin = T)$Odata
psi0.Y <- mean(indPop.iid.cA.bY$Y)  # 0.291398
psi0.Ygstar <- mean(indPop.iid.cA.bY$Y.gstar)  # 0.316274
indSample.iid.cA.bY <- get.iid.dat.Acont(ndata = ndata, rndseed = rndseed, truncBD = truncBD, shift.val = shift.val, is.Y.bin = T)$Odata
indSample.iid.cA.bY <- indSample.iid.cA.bY[, c("W1", "W2", "W3", "W4", "A", "Y", "Y.gstar", "trunc.A.gstar")]
indSample.iid.cA.bY_list <- list(indSample.iid.cA.bY = indSample.iid.cA.bY, psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar,
                                 truncBD = truncBD, shift.val = shift.val, rndseed = rndseed)
save(indSample.iid.cA.bY_list, file="indSample.iid.cA.bY_list.Rda")


#### Data 2. One continuous A with Continuous Y
shift.val <- 2
indPop.iid.cA.cY <- get.iid.dat.Acont(ndata = 1000000, rndseed = rndseed, truncBD = truncBD, shift.val = shift.val, is.Y.bin = F)$Odata
psi0.Y <- mean(indPop.iid.cA.cY$Y)  # 3.309084
psi0.Ygstar <- mean(indPop.iid.cA.cY$Y.gstar)  # 3.50856
indSample.iid.cA.cY <- get.iid.dat.Acont(ndata = ndata, rndseed = rndseed, truncBD = truncBD, shift.val = shift.val, is.Y.bin = F)$Odata
indSample.iid.cA.cY <- indSample.iid.cA.cY[, c("W1", "W2", "W3", "W4", "A", "Y", "Y.gstar", "trunc.A.gstar")]
indSample.iid.cA.cY_list <- list(indSample.iid.cA.cY = indSample.iid.cA.cY, psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar,
                                 truncBD = truncBD, shift.val = shift.val, rndseed = rndseed)
                                 # call = "get.iid.dat.Acont(ndata=10000, rndseed=12345, truncBD=10, shift.val=2, is.Y.bin=F)$Odata")
save(indSample.iid.cA.cY_list, file="indSample.iid.cA.cY_list.Rda")
