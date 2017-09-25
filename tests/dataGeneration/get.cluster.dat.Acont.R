# -------------------------------------------------------------------------------------------------------
# DATA SET 4. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN CLUSTER DATA (IID within Cluster)
# -------------------------------------------------------------------------------------------------------
# Fitting continuous exposure by logistic regression (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (E1, E2, W1, W2);
# -------------------------------------------------------------------------------------------------------

get.cluster.dat.Acont <- function(id, n.incluster = 10000, rndseed = NULL, is.Y.bin = T, truncBD = 10, shift.val = 2, 
                                  working.model = T, timevarying = F) {
  set.seed(rndseed)
  E1 <- runif(n = 1, min = 0, max = 1)
  E2 <- sample(x = c(0, 0.5, 1), size = 1)
  W1 <- rbinom(n = n.incluster, size = 1, prob = plogis(- 0.7 + 1.2 * E1 - 2 * E2))
  UW2 <- MASS::mvrnorm(n = n.incluster, mu = rep(0, 2), Sigma = matrix(c(1, 0.7, 0.7, 1), ncol = 2))
  W2 <- 0.8 * E1 + 0.6 * E2 + 1 * UW2[, 1] # + 1 * UW2[, 2]
  if (timevarying) {
    A <- rbinom(n = 1, size = 1, prob = plogis(- 2.8 + 3.5 * E1 + 0.4 * E2 + 6 * mean(W1)))  # community-level A
  }
  if (working.model) {  # working.model holds
    if (is.Y.bin) {
      Y <- rbinom(n = n.incluster, size = 1, prob = plogis(- 2 + 2 * A + 2 * E1 + 0.6 * E2 + W1 + 0.6 * W2))
      Y1 <- rbinom(n = n.incluster, size = 1, prob = plogis(- 2 + 2 * 1 + 2 * E1 + 0.6 * E2 + W1 + 0.6 * W2))
      Y0 <- rbinom(n = n.incluster, size = 1, prob = plogis(- 2 + 2 * 0 + 2 * E1 + 0.6 * E2 + W1 + 0.6 * W2))
    } else {
      # finishing ....
    }
  } else {  # working.model fails
    if (is.Y.bin) {
      Y <- rbinom(n = n.incluster, size = 1, prob = plogis(-3 + 1 * A + 0.2 * E1 + 0.2 * E2 + 0.1 * W1 + 0.2 * W2 + 3 * mean(W1) + 1.5 * mean(W2)))
      Y1 <- rbinom(n = n.incluster, size = 1, prob = plogis(-3 + 1 * 1 + 0.2 * E1 + 0.2 * E2 + 0.1 * W1 + 0.2 * W2 + 3 * mean(W1) + 1.5 * mean(W2)))
      Y0 <- rbinom(n = n.incluster, size = 1, prob = plogis(-3 + 1 * 0 + 0.2 * E1 + 0.2 * E2 + 0.1 * W1 + 0.2 * W2 + 3 * mean(W1) + 1.5 * mean(W2)))
    } else {
      # finishing ....
    }
  }
  return(data.frame(cbind(id = id, E1 = E1, E2 = E2, W1 = W1, W2 = W2, A = A, Y = Y, Y1 = Y1, Y0 = Y0)))
}
