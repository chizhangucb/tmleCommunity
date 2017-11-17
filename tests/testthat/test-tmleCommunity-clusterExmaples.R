context("Test tmleCommunity on community-level exposure")

# ---------------------------------------------------------------------------------
# TEST SET 1. TESTS FOR FITTING BINARY EXPOSURE A IN HIERARCHICAL DATA
# ---------------------------------------------------------------------------------
get.cluster.dat.Abin <- function(id, n.ind = 10000, rndseed = NULL, is.Y.bin = TRUE) {
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
  # working.model fails since other individuals' W affect other's Y
  Y <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.8 * A + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 
                                                 - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
  Y1 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.8 * 1 + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 
                                                  - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
  Y0 <- rbinom(n = n.ind, size = 1, prob = plogis(- 0.2 + 0.8 * 0 + 0.5 * E1 + 0.2 * E2 + 0.3 * W1 
                                                  - 0.15 * W3 + 0.75 * mean(W2) + 0.4 * mean(W3)))
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, Y, Y1, Y0)))
}

`%+%` <- function(a, b) paste0(a, b)
J <- 1000
n.ind <- 50
rndseed <- 12345
comSample.wmF.bA.bY <- get.fullDat.Abin(J = J, n.ind = n.ind, rndseed = rndseed, is.Y.bin = TRUE, working.model = FALSE)
comSample.wmF.bA.bY <- comSample.wmF.bA.bY[, c("id", "E1", "E2", "W1", "W2", "W3", "A", "Y")]
N <- NROW(comSample.wmF.bA.bY)
Qform.corr <- "Y ~ E1 + E2 + W1 + W2 + W3 + A" # correct Q form
gform.corr <- "A ~ E1 + E2 + W1"  # correct g form

#*************************************** 
## Test 1.1 speed.glm & glm 
#*************************************** 
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)

## Test 1.1.1 Community-level analysis without a pooled individual-level regression on outcome
test_that("fit TMLE for bin, community-level A, with community-level analysis without a pooled Q", {
   tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                               WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                               f_gstar2 = 0L, community.step = "community_level", 
                               communityID = "id", pooled.Q = FALSE, Qform = Qform.corr, 
                               hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 
  expect_equal(estimates["tmle", ],  , tolerance = 1e-6))
  expect_equal(estimates["iptw", ],  , tolerance = 1e-6))
  expect_equal(estimates["gcomp", ],  , tolerance = 1e-6))
  expect_that()
})

## Test 1.1.2 Community-level analysis with a pooled individual-level regression on outcome
test_that("fit TMLE for bin, community-level A, with community-level analysis with a pooled Q", {
   tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                               WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                               f_gstar2 = 0L, community.step = "community_level", 
                               communityID = "id", pooled.Q = TRUE, Qform = Qform.corr, 
                               hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 
  expect_equal(estimates["tmle", ],  , tolerance = 1e-6))
  expect_equal(estimates["iptw", ],  , tolerance = 1e-6))
  expect_equal(estimates["gcomp", ],  , tolerance = 1e-6))
  expect_that()
})

## Test 1.1.3 Individual-level analysis with both individual-level outcome and treatment mechanisms
test_that("fit TMLE for bin, community-level A, with individual-level analysis", {
   tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                               WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                               f_gstar2 = 0L, community.step = "individual_level", 
                               communityID = "id", pooled.Q = TRUE, Qform = Qform.corr, 
                               hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 
  expect_equal(estimates["tmle", ],  , tolerance = 1e-6))
  expect_equal(estimates["iptw", ],  , tolerance = 1e-6))
  expect_equal(estimates["gcomp", ],  , tolerance = 1e-6))
  expect_that()
})

#***************************************************************************************

# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN HIERARCHICAL DATA
# ---------------------------------------------------------------------------------
