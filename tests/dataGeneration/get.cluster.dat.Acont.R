# -------------------------------------------------------------------------------------------------------
# DATA SET 4. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN CLUSTER DATA (IID within Cluster)
# -------------------------------------------------------------------------------------------------------
# Fitting continuous exposure by logistic regression (One A with Binary/ Continuous Y)
# individual exposure is normal with mu for each observation being a function of (E1, E2, W1, W2);
# -------------------------------------------------------------------------------------------------------

get.cluster.dat.Acont <- function(id, n.incluster = 1000, rndseed = NULL, is.Y.bin = T, truncBD = 5, shift.val = 1, 
                                  working.model = T, timevarying = F) {
  set.seed(rndseed)
  # Construct community-level & individual-level baseline covariates E, W 
  E1 <- runif(n = 1, min = 0, max = 1)
  E2 <- sample(x = c(0, 0.25, 0.75, 1), size = 1)
  W1 <- rbinom(n = n.incluster, size = 1, prob = plogis(- 0.4 + 1.2 * E1 - 1.3 * E2))
  UW2 <- MASS::mvrnorm(n = n.incluster, mu = rep(0, 2), Sigma = matrix(c(1, 0.6, 0.6, 1), ncol = 2))
  W2 <- 0.8 * E1 - 0.4 * E2 + 1 * UW2[, 1] 
  W3 <- rbinom(n = n.incluster, size = 1, prob = 0.5)
  
  # Construct continuous exposure; time-varying within community if timevarying = T
  if (timevarying) { # A is time-varying within community 
    A.mu <- - 1.2 + 0.8 * E1 + 0.21 * E2 + 3 * W1 - 0.7 * W2 + 1.3 * W3
    A <- rnorm(n = n.incluster, mean = A.mu, sd = 1)
  } else {  # A is constant within community 
    A.mu <- - 1.2 + 0.8 * E1 + 0.21 * E2 + 3 * mean(W1) - 0.7 * mean(W2) + 1.3 * mean(W3)
    A <- rnorm(n = 1, mean = A.mu, sd = 1)
  }
  untrunc.A.gstar <- A + shift.val
  r.new.A <- exp(2.8 * shift.val * (untrunc.A.gstar - A.mu - shift.val / 5))
  trunc.A.gstar <- ifelse(r.new.A > truncBD, A, untrunc.A.gstar)
  
  # Construct outcome
  if (working.model) {  # when working.model holds
    if (is.Y.bin) {
      Y <- rbinom(n = n.incluster, size = 1, prob = plogis(- 1.7 + 2.4 * A + 0.5 * E1 - 1.2 * E2 + 0.8 * W1 + 1 * W2 - 0.4 * W3))
      Y.gstar <- rbinom(n = n.incluster, size = 1, prob = plogis(- 1.7 + 2.4 * trunc.A.gstar + 0.5 * E1 - 1.2 * E2 + 0.8 * W1 + 1 * W2 - 0.4 * W3))
    } else {
      # finishing ....
    }
  } else {  # when working.model fails
    if (is.Y.bin) {
      Y <- rbinom(n = n.incluster, size = 1, 
                  prob = plogis(- 2 + 1 * A + 0.3 * E1 + 1.1 * E2 - 0.4 * W1 + 0.1 * W2 + 1.8 * mean(W1) + 2 * mean(W2)) + 3 * mean(W3))
      Y.gstar <- rbinom(n = n.incluster, size = 1, prob = plogis(- 2 + 1 * trunc.A.gstar + 0.3 * E1 + 1.1 * E2 - 0.4 * W1 +
                                                                   0.1 * W2 + 1.8 * mean(W1) + 2 * mean(W2)) + 3 * mean(W3))
    } else {
      # finishing ....
    }
  }
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, trunc.A.gstar, Y, Y.gstar)))
}

