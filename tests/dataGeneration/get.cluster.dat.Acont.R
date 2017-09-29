# -------------------------------------------------------------------------------------------------------
# DATA SET 3. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN CLUSTER DATA (IID within Cluster)
# -------------------------------------------------------------------------------------------------------
# Fitting continuous exposure by logistic regression (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (E1, E2, W1, W2);
# -------------------------------------------------------------------------------------------------------

get.cluster.dat.Acont <- function(id, n.ind = 1000, rndseed = NULL, is.Y.bin = T, truncBD = 5, shift.val = 1, 
                                  working.model = T, timevarying = T) {
  set.seed(rndseed)
  # Construct community-level & individual-level baseline covariates E, W 
  E1 <- runif(n = 1, min = 0, max = 1)
  E2 <- sample(x = c(0.2, 0.4, 0.6, 0.8), size = 1)
  W1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.4 + 1.2 * E1 - 1.3 * E2))
  UW2 <- MASS::mvrnorm(n = n.ind, mu = rep(0, 2), Sigma = matrix(c(1, 0.6, 0.6, 1), ncol = 2))
  W2 <- 0.8 * E1 - 0.4 * E2 + 1 * UW2[, 1] 
  W3 <- sample(x = 1:12, size = n.ind, replace = TRUE) / 10
  
  # Construct continuous exposure; time-varying within community if timevarying = T
  if (timevarying) { # A is time-varying within community 
    A.mu <- - 1.2 + 0.8 * E1 + 0.21 * E2 + 3 * W1 - 0.7 * W2 + 0.3 * W3
    A <- rnorm(n = n.ind, mean = A.mu, sd = 1)
  } else {  # A is constant within community 
    A.mu <- - 1.2 + 0.8 * E1 + 0.21 * E2 + 3 * mean(W1) - 0.7 * mean(W2) + 0.3 * mean(W3)
    A <- rnorm(n = 1, mean = A.mu, sd = 1)
  }
  untrunc.A.gstar <- A + shift.val
  r.new.A <- exp(1.8 * shift.val * (untrunc.A.gstar - A.mu - shift.val / 4))
  trunc.A.gstar <- ifelse(r.new.A > truncBD, A, untrunc.A.gstar)
  
  # Construct outcome
  if (working.model) {  # when working.model holds
    if (is.Y.bin) {
      Y <- rbinom(n = n.ind, size = 1, prob = plogis(- 1.7 + 3.4 * A + 0.5 * E1 - 1.2 * E2 + 0.8 * W1 + 1 * W2 - 0.1 * W3))
      Y.gstar <- rbinom(n = n.ind, size = 1, prob = plogis(- 1.7 + 3.4 * trunc.A.gstar + 0.5 * E1 - 1.2 * E2 + 0.8 * W1 + 1 * W2 - 0.1 * W3))
    } else {
      # finishing ....
    }
  } else {  # when working.model fails
    if (is.Y.bin) {
      Y <- rbinom(n = n.ind, size = 1, prob = plogis(- 5 + 1 * A - 0.2 * E1 + 1.1 * E2 - 0.4 * W1 + 0.2 * W2 - 0.3 * W3 + 
                                                       5.8 * mean(W1) + 3.1 * mean(W2) - 1 * mean(W3)))
      Y.gstar <- rbinom(n = n.ind, size = 1, prob = plogis(- 5 + 1 * trunc.A.gstar - 0.2 * E1 + 1.1 * E2 - 0.4 * W1 + 0.2 * W2 - 
                                                             0.3 * W3 + 5.8 * mean(W1) + 3.1 * mean(W2) - 1 * mean(W3)))
    } else {
      # finishing ....
    }
  }
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, trunc.A.gstar, Y, Y.gstar)))
}


get.fullDat.Acont <- function(J, n.ind, rndseed = NULL, is.Y.bin = TRUE, truncBD = 5, shift.val = 1, working.model = TRUE, 
                              timevarying = TRUE, n.ind.fix = FALSE, onlyYkeep = FALSE, verbose = TRUE) {
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
    cluster.data.j <- get.cluster.dat.Acont(id = j, n.ind = n.ind[j], rndseed = rndseed, is.Y.bin = is.Y.bin, truncBD = truncBD, 
                                            shift.val = shift.val, working.model = working.model, timevarying = timevarying) 
    if (onlyYkeep) {
      Y <- c(Y, cluster.data.j[, "Y"]); Y.gstar <- c(Y.gstar, cluster.data.j[, "Y.gstar"]) 
    } else {
      full.data <- rbind(full.data, cluster.data.j)
    }
  }  
  ifelse(onlyYkeep, return(data.frame(cbind(Y, Y.gstar))), return(full.data))
}

J <- 1000
n.ind <- 100
rndseed <- 12345
truncBD <- 5
shift.val <- 1

#### Data 1. One A with Binary Y, when working model fails & A is time-varying
popDat.wmF.tvcontA.binY <- get.fullDat.Acont(J = 4000, n.ind = 4000, rndseed = NULL, truncBD = truncBD, shift.val = shift.val, 
                                             is.Y.bin = TRUE, working.model = FALSE, timevarying = TRUE, n.ind.fix = FALSE, onlyYkeep = T)
psi0.Y <- mean(popDat.wmF.tvcontA.binY$Y)  # 0.2406796
psi0.Ygstar <- mean(popDat.wmF.tvcontA.binY$Y.gstar)  # 0.2964227
dat_com.wmF.tvcontA.binY <- get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = rndseed, truncBD = truncBD, shift.val = shift.val, 
                                              is.Y.bin = TRUE, working.model = FALSE, timevarying = TRUE, n.ind.fix = FALSE)
sampleDat_com.wmF.tvcontA.binY <- list(dat_com.wmF.tvcontA.binY = dat_com.wmF.tvcontA.binY, psi0.Y = psi0.Y, psi0.Ygstar = psi0.Ygstar,
                                       truncBD = truncBD, shift.val = shift.val, rndseed = rndseed)
save(sampleDat_com.wmF.tvcontA.binY, file = "sampleDat_com.wmF.tvcontA.binY.Rda")

#### Data 2. One A with Binary Y, when working model holds & A is time-varying

