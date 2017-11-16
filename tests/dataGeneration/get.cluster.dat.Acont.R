# -------------------------------------------------------------------------------------------------------
# DATA SET 3. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN CLUSTER DATA (IID within Cluster)
# -------------------------------------------------------------------------------------------------------
# Fitting continuous exposure by logistic regression (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (E1, E2, W1, W2);
# -------------------------------------------------------------------------------------------------------

get.cluster.dat.Acont <- function(id, n.ind = 1000, is.Y.bin = TRUE, truncBD = 5, shift.val = 1, 
                                  rndseed = NULL, working.model = TRUE) {
  set.seed(rndseed)
  # Construct community-level & individual-level baseline covariates E, W 
  E1 <- runif(n = 1, min = 0, max = 1)
  E2 <- sample(x = c(0.2, 0.4, 0.6, 0.8), size = 1)
  W1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.4 + 1.2 * E1 - 1.3 * E2))
  W2_mean <- 1 - 0.8 * E1 - 0.4 * E2 
  W3_mean <- 0.5 + 0.2 * E1
  W2W3 <- MASS::mvrnorm(n = n.ind, mu = c(W2_mean, W3_mean), Sigma = matrix(c(1, 0.6, 0.6, 1), ncol = 2))
  W2 <- W2W3[, 1]
  W3 <- W2W3[, 2] 
  A.mu <- - 1.2 + 0.8 * E1 + 0.21 * E2 + 3 * mean(W1) - 0.7 * mean(W2) + 0.3 * mean(W3)
  A <- rnorm(n = 1, mean = A.mu, sd = 1)
  untrunc.A.gstar <- A + shift.val
  r.new.A <- exp(1.5 * shift.val * (untrunc.A.gstar - A.mu - shift.val / 4))
  trunc.A.gstar <- ifelse(r.new.A > truncBD, A, untrunc.A.gstar)
  
  # Construct outcome
  if (working.model) {  # when working.model holds
    if (is.Y.bin) {  # Create binary Y
      Y <- rbinom(n = n.ind, size = 1, 
                  prob = plogis(- 1.7 + 1.2 * A + 0.5 * E1 - 1.2 * E2 + 0.7 * W1 + 1.3 * W2 - 0.4 * W3))
      Y.gstar <- rbinom(n = n.ind, size = 1,
                        prob = plogis(- 1.7 + 1.2 * trunc.A.gstar + 0.5 * E1 - 1.2 * E2 + 0.7 * W1 + 1.3 * W2 - 0.4 * W3))
    } else {  # Create continuous Y
      Y <- rnorm(n = 1, mean =  -2 + 2 * A + 2 * E1 + 0.6 * E2 + W1 + 0.6 * W2 - 2.1 * W3, sd = 1)
      Y <- rnorm(n = 1, mean =  -2 + 2 * trunc.A.gstar + 2 * E1 + 0.6 * E2 + W1 + 0.6 * W2 - 2.1 * W3, sd = 1)
    }
  } else {  # when working.model fails
    if (is.Y.bin) {  # Create binary Y
      Y <- rbinom(n = n.ind, size = 1, 
                  prob = plogis(- 5 + 1 * A - 0.2 * E1 + 1.1 * E2 - 0.4 * W1 + 0.2 * W2 - 0.4 * W3 + 
                                  5.8 * mean(W1) + 3.1 * mean(W2) - 1 * mean(W3)))
      Y.gstar <- rbinom(n = n.ind, size = 1, 
                        prob = plogis(- 5 + 1 * trunc.A.gstar - 0.2 * E1 + 1.1 * E2 - 0.4 * W1 + 0.2 * W2 - 
                                        0.4 * W3 + 5.8 * mean(W1) + 3.1 * mean(W2) - 1 * mean(W3)))
    } else {  # Create continuous Y
      # finishing ....
    }
  }
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, trunc.A.gstar, Y, Y.gstar)))
}

get.fullDat.Acont <- function(J, n.ind, rndseed = NULL, is.Y.bin = TRUE, truncBD = 5, shift.val = 1, 
                              working.model = TRUE, n.ind.fix = FALSE, onlyYkeep = FALSE, verbose = TRUE) {
  set.seed(rndseed)
  if (n.ind.fix) {
    n.ind <- rep(n.ind, J)
  } else {  # don't fix the number of individuals in each community
    n.ind <- round(rnorm(J, n.ind, 10))  
    n.ind[n.ind <= 0] <- n.ind  # set to n.ind if any generated number less than 1
  }
  
  if (onlyYkeep) {
    message("Only the unshifted & shifted outcomes are useful, only just keep these two")
    Y <- Y.gstar <- NULL
  } else {
    full.data <- NULL
  }
  
  for(j in 1:J) {
    if (verbose) message("#### generating " %+% j %+% "th cluster ####")
    cluster.data.j <- get.cluster.dat.Acont(id = j, n.ind = n.ind[j], rndseed = rndseed, is.Y.bin = is.Y.bin, 
                                            truncBD = truncBD, shift.val = shift.val, working.model = working.model) 
    if (onlyYkeep) {
      Y <- c(Y, cluster.data.j[, "Y"]); Y.gstar <- c(Y.gstar, cluster.data.j[, "Y.gstar"]) 
    } else {
      full.data <- rbind(full.data, cluster.data.j)
    }
  }  
  ifelse(onlyYkeep, return(data.frame(cbind(Y, Y.gstar))), return(full.data))
}

J <- 1000
n.ind <- 50
rndseed <- 12345
truncBD <- 5
shift.val <- 1

#### Data 1. One continuous, community-level A with binary Y, when working model holds
comPop.wmT.cA.bY <- get.fullDat.Acont(J = 4000, n.ind = 4000, rndseed = NULL, is.Y.bin = TRUE, truncBD = truncBD,
                                      shift.val = shift.val, working.model = TRUE, onlyYkeep = TRUE)
mean(comPop.wmT.bA.bY$Y1) # 0.5150197
mean(comPop.wmT.bA.bY$Y0) # 0.4113038
psi0.Y <- mean(comPop.wmT.cA.bY$Y)  # 0.2472083
psi0.Ygstar <- mean(comPop.wmT.cA.bY$Y.gstar)  # 0.3123667

comSample.wmT.cA.bY <- get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = rndseed, truncBD = truncBD, 
                                         shift.val = shift.val, is.Y.bin = TRUE, working.model = TRUE)
comSample.wmT.cA.bY <- comSample.wmT.cA.bY[, c("id", "E1", "E2", "W1", "W2", "W3", "A", "Y")]
comSample.wmT.cA.bY_list <- list(comSample.wmT.cA.bY = comSample.wmT.cA.bY, psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar,
                                 truncBD = truncBD, shift.val = shift.val, rndseed = rndseed)
save(comSample.wmT.cA.bY_list, file = "comSample.wmT.cA.bY_list.Rda")

#### Data 2. One continuous, community-level A with binary Y, when working model fails
comPop.wmT.cA.bY <- get.fullDat.Acont(J = 4000, n.ind = 4000, rndseed = NULL, is.Y.bin = TRUE, truncBD = truncBD,
                                      shift.val = shift.val, working.model = FALSE, onlyYkeep = TRUE)
psi0.Y <- mean(comPop.wmT.cA.bY$Y)  # 0.2860212
psi0.Ygstar <- mean(comPop.wmT.cA.bY$Y.gstar)  # 0.3340003

comSample.wmF.cA.bY <- get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = rndseed, truncBD = truncBD, 
                                         shift.val = shift.val, is.Y.bin = TRUE, working.model = FALSE)
comSample.wmF.cA.bY <- comSample.wmF.cA.bY[, c("id", "E1", "E2", "W1", "W2", "W3", "A", "Y")]
comSample.wmF.cA.bY_list <- list(comSample.wmF.cA.bY = comSample.wmF.cA.bY, psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar,
                                 truncBD = truncBD, shift.val = shift.val, rndseed = rndseed)
save(comSample.wmF.cA.bY_list, file = "comSample.wmF.cA.bY_list.Rda")
