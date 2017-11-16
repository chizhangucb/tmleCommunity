# -------------------------------------------------------------------------------------------------------
# DATA SET 4. TESTS FOR FITTING COMMUNITY-LEVEL BINARY EXPOSURE A IN CLUSTER DATA (IID within Cluster)
# -------------------------------------------------------------------------------------------------------
# Fitting binary exposure by logistic regression (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (E1, E2, W1, W2);
# -------------------------------------------------------------------------------------------------------

get.cluster.dat.Abin <- function(id, n.ind = 10000, rndseed = NULL, is.Y.bin = TRUE, working.model = TRUE) {
  set.seed(rndseed)
  E1 <- runif(n = 1, min = 0, max = 1)
  E2 <- sample(x = c(0, 0.2, 0.4, 0.8, 1), size = 1)
  W1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.13 + 1.5 * E1 - 0.6 * E2))
  W2_mean <- - 0.6 + 0.8 * E1 - 0.4 * E2 
  W3_mean <- 0.5 + 0.2 * E1
  W2W3 <- MASS::mvrnorm(n = n.ind, mu = c(W2_mean, W3_mean), Sigma = matrix(c(1, 0.6, 0.6, 1), ncol = 2))
  W2 <- W2W3[, 1]
  W3 <- W2W3[, 2] 
  A <- rbinom(n = 1, size = 1, prob = plogis(- 2.4 + 3 * E1 + 0.4 * E2 + 3 * mean(W1)))  # community-level A
  if (working.model) {  # working.model holds
    if (is.Y.bin) {  # Y binary
      Y <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.5 * A + 0.4 * E1 + 0.2 * E2 + W2 - 0.2 * W3))
      Y1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.5 * 1 + 0.4 * E1 + 0.2 * E2 + W2 - 0.2 * W3))
      Y0 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.55 * 0 + 0.4 * E1 + 0.2 * E2 + W2 - 0.2 * W3))
    } else {  # Y continuous 
      Y <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.4 * W2 + 0.2 * W3, sd = 1)
      Y1 <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.4 * W2 + 0.2 * W3, sd = 1)
      Y0 <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.4 * W2 + 0.2 * W3, sd = 1)
    }
  } else {  # working.model fails
    if (is.Y.bin) {  # Y binary
      Y <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.3 * A + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
      Y1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.3 * 1 + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
      Y0 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.3 * 0 + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
    } else {  # Y continuous 
      Y <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3), sd = 1)
      Y1 <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3), sd = 1)
      Y0 <- rnorm(n = n.ind, mean = 1 + 0.3 * A - 0.2 * E1 + 0.4 * E2 + 0.3 * W1 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3), sd = 1)
    }
  }
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, Y, Y1, Y0)))
}

get.fullDat.Abin <- function(J, n.ind, rndseed = NULL, is.Y.bin = TRUE, working.model = TRUE, 
                             n.ind.fix = FALSE, onlyYkeep = FALSE, verbose = TRUE) {
  set.seed(rndseed)
  if (n.ind.fix) {
    n.ind <- rep(n.ind, J)
  } else {  # don't fix the number of individuals in each community
    n.ind <- round(rnorm(J, n.ind, 10))  
    n.ind[n.ind <= 0] <- n.ind  # set to n.ind if any generated number less than 1
  }
  
  if (onlyYkeep) {
    message("Only the observed & (two) post-intervened outcomes are useful, only just keep these three")
    Y <- Y1 <- Y0 <- NULL
  } else {
    full.data <- NULL
  }
  
  for(j in 1:J) {
    if (verbose) message("#### generating " %+% j %+% "th cluster ####")
    cluster.data.j <- get.cluster.dat.Abin(id = j, n.ind = n.ind[j], is.Y.bin = is.Y.bin, working.model = working.model) 
    if (onlyYkeep) {
      Y <- c(Y, cluster.data.j[, "Y"])
      Y1 <- c(Y1, cluster.data.j[, "Y1"])
      Y0 <- c(Y0, cluster.data.j[, "Y0"])
    } else {
      full.data <- rbind(full.data, cluster.data.j)
    }
  }  
  full.data$id <- as.integer(full.data$id)
  ifelse(onlyYkeep, return(data.frame(cbind(Y, Y1, Y0))), return(full.data))
}

J <- 1000
n.ind <- 50
rndseed <- 12345

#### Data 1. One binary, community-level A with Binary Y, when working model holds
comPop.wmT.bA.bY <- get.fullDat.Abin(J = 4000, n.ind = 4000, rndseed = NULL, is.Y.bin = T, working.model = T, n.ind.fix = F, onlyYkeep = T)
mean(comPop.wmT.bA.bY$Y1) # 0.5150197
mean(comPop.wmT.bA.bY$Y0) # 0.4113038
psi0.Y <- mean(comPop.wmT.bA.bY$Y1) - mean(comPop.wmT.bA.bY$Y0) # 0.103716

comSample.wmT.bA.bY <- get.fullDat.Abin(J = J, n.ind = n.ind, rndseed = rndseed, is.Y.bin = T, working.model = T, n.ind.fix = F)
comSample.wmT.bA.bY <- comSample.wmT.bA.bY[, c("id", "E1", "E2", "W1", "W2", "W3", "A", "Y")]
comSample.wmT.bA.bY_list <- list(comSample.wmT.bA.bY = comSample.wmT.bA.bY, psi0.Y = psi0.Y, rndseed = rndseed)
comSample.wmT.bA.bY_list$psi0.Y  # 0.103716
save(comSample.wmT.bA.bY_list, file = "comSample.wmT.bA.bY_list.rda")

#### Data 2. One binary, community-level A with Binary Y, when working model fails
comPop.wmF.bA.bY <- get.fullDat.Abin(J = 4000, n.ind = 4000, rndseed = NULL, is.Y.bin = T, working.model = F, n.ind.fix = F, onlyYkeep = T)
mean(comPop.wmF.bA.bY$Y1) # 
mean(comPop.wmF.bA.bY$Y0) # 
psi0.Y <- mean(comPop.wmF.bA.bY$Y1) - mean(comPop.wmF.bA.bY$Y0) # 

comSample.wmF.bA.bY <- get.fullDat.Abin(J = J, n.ind = n.ind, rndseed = rndseed, is.Y.bin = T, working.model = T, n.ind.fix = F)
comSample.wmF.bA.bY_list <- list(comSample.wmF.bA.bY = comSample.wmF.bA.bY, psi0.Y = psi0.Y, rndseed = rndseed)
save(comSample.wmF.bA.bY_list, file = "comSample.wmF.bA.bY_list.rda")
