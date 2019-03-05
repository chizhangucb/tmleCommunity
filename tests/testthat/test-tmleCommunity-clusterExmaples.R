context("Test tmleCommunity on community-level exposure")

# ---------------------------------------------------------------------------------
# TEST SET 1. TESTS FOR FITTING BINARY EXPOSURE A IN HIERARCHICAL DATA
# ---------------------------------------------------------------------------------
get.cluster.dat.Abin <- function(id, n.ind = 10000, rndseed = NULL) {
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

get.fullDat.Abin <- function(J, n.ind, rndseed = NULL, n.ind.fix = FALSE, verbose = TRUE) {
  set.seed(rndseed)
  if (n.ind.fix) {
    n.ind <- rep(n.ind, J)
  } else {  # don't fix the number of individuals in each community
    n.ind <- round(rnorm(J, n.ind, 10))  
    n.ind[n.ind <= 0] <- n.ind  # set to n.ind if any generated number less than 1
  }
  full.data <- NULL

  for(j in 1:J) {
    if (verbose) message("#### generating " %+% j %+% "th cluster ####")
    cluster.data.j <- get.cluster.dat.Abin(id = j, n.ind = n.ind[j], rndseed = eval(rndseed + j)) 
    full.data <- rbind(full.data, cluster.data.j)
  }  
  full.data$id <- as.integer(full.data$id)
  return(full.data)
}

`%+%` <- function(a, b) paste0(a, b)
J <- 1000
n.ind <- 50
rndseed <- 12345
comSample.wmF.bA.bY <- get.fullDat.Abin(J = J, n.ind = n.ind, rndseed = rndseed)
mean(comSample.wmF.bA.bY$Y1) - mean(comSample.wmF.bA.bY$Y0) # 0.1753572
comSample.wmF.bA.bY <- comSample.wmF.bA.bY[, c("id", "E1", "E2", "W1", "W2", "W3", "A", "Y")]
N <- NROW(comSample.wmF.bA.bY)
Qform.corr <- "Y ~ E1 + E2 + W1 + W2 + W3 + A" # correct Q form
gform.corr <- "A ~ E1 + E2 + W1"  # correct g form

#*************************************** 
## Test 1.1 Different analysis methods
#*************************************** 
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)

## Test 1.1.1 Community-level analysis without a pooled individual-level regression on outcome
test_that("fit TMLE for bin, community-level A, with community-level analysis without a pooled Q", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.step = "community_level", 
                              communityID = "id", pooled.Q = FALSE)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  expect_equal(estimates["tmle", ], 0.1773423, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1899487, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1813821, tolerance = 1e-6)
})

## Test 1.1.2 Community-level analysis with a pooled individual-level regression on outcome
test_that("fit TMLE for bin, community-level A, with community-level analysis with a pooled Q", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.step = "community_level", 
                              communityID = "id", pooled.Q = TRUE)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  expect_equal(estimates["tmle", ], 0.1784616, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1957239, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1838956, tolerance = 1e-6)
})

# Test 1.1.3 Individual-level analysis with both individual-level outcome and treatment mechanisms
# the bias of estimation is larger than community-level analysis since the working model fails
test_that("fit TMLE for bin, community-level A, with individual-level analysis", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.step = "individual_level", communityID = "id")
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  variances <- tmleCom_res$ATE$vars
  expect_equal(estimates["tmle", ], 0.1835877, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1854157, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1825362, tolerance = 1e-6)
  expect_equal(as.vector(variances[, 1]), c(3.899418e-05, 2.864193e-03, 3.899418e-05))
})

# Test 1.1.4 Individual-level analysis that assumes no hierarchical structure
test_that("Treat hierarchical data as non-hierarchical and use usual TMLE", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.step = "NoCommunity")
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  variances <- tmleCom_res$ATE$vars
  expect_equal(estimates["tmle", ], 0.1835877, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1854157, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1825362, tolerance = 1e-6)
  expect_equal(as.vector(variances[, 1]), c(3.217030e-05, 8.865422e-05, 3.217030e-05))
  
  expect_message(
    tmleCom_res2 <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                                 WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                                 f_gstar2 = 0L, community.step = "individual_level", 
                                 communityID = NULL),
    regexp = "Lack of 'communityID' forces the algorithm to automatically " %+%
      "pool data over all communities and treat it as non-hierarchical dataset"
  )
  estimates2 <- tmleCom_res2$ATE$estimates  # psi0 = 0.1768159
  expect_equal(estimates, estimates2)
})

## Test 1.1.5 Unweighted covariate-based TMLE, using community-level analysis without a pooled Q
test_that("fit unweighted covariate-based TMLE for bin, community-level A", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.step = "community_level", 
                              communityID = "id", TMLE.targetStep = "tmle.covariate")
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  expect_equal(estimates["tmle", ], 0.1764060, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1899487, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1813821, tolerance = 1e-6)
})

## Test 1.1.6 Weigh each community equally, using community-level analysis without a pooled Q
test_that("fit TMLE for bin, community-level A, while weighing each community equally", {
  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
                              f_gstar2 = 0L, community.wts = "equal.community",
                              community.step = "community_level", communityID = "id")
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
  expect_equal(estimates["tmle", ], 0.1766701, tolerance = 1e-6)
  expect_equal(estimates["iptw", ], 0.1907880, tolerance = 1e-6)
  expect_equal(estimates["gcomp", ], 0.1809624, tolerance = 1e-6)
})

# Test 1.1.7 Stratified TMLE that fits a separate outcome (exposure) mechanism for each community
#test_that("fit stratified TMLE for bin, community-level A", {
#  tmleCom_Options(Qestimator = "glm__glm", gestimator = "glm__glm", maxNperBin = N)
#  tmleCom_res <-tmleCommunity(data = comSample.wmF.bA.bY, Ynode = "Y", Anodes = "A",
#                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, 
#                              f_gstar2 = 0L, community.step = "perCommunity", communityID = "id")
#  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.1768159
#  expect_equal(estimates["tmle", ], 0.03036185, tolerance = 1e-6)
#  expect_equal(estimates["gcomp", ], 0.00000000, tolerance = 1e-6)
#  expect_equal(as.vector(variances[, 1]), c(2.08973e-11, 3.41511e-04, 2.08973e-11))
#})

# ---------------------------------------------------------------------------------
# TEST SET 2. TESTS FOR FITTING CONTINUOUS EXPOSURE A IN HIERARCHICAL DATA
# ---------------------------------------------------------------------------------
get.cluster.dat.Acont <- function(id, n.ind = 1000, truncBD = 5, shift.val = 1, rndseed = NULL) {
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
  Y <- rbinom(n = n.ind, size = 1, 
              prob = plogis(- 1.7 + 1.2 * A + 0.5 * E1 - 1.2 * E2 + 0.7 * W1 + 1.3 * W2 - 0.4 * W3))
  Y.gstar <- rbinom(n = n.ind, size = 1,
                    prob = plogis(- 1.7 + 1.2 * trunc.A.gstar + 0.5 * E1 - 1.2 * E2 + 0.7 * W1 + 1.3 * W2 - 0.4 * W3))
  return(data.frame(cbind(id, E1, E2, W1, W2, W3, A, trunc.A.gstar, Y, Y.gstar)))
}

get.fullDat.Acont <- function(J, n.ind, rndseed = NULL, truncBD = 5, shift.val = 1, 
                              n.ind.fix = FALSE, verbose = TRUE) {
  set.seed(rndseed)
  if (n.ind.fix) {
    n.ind <- rep(n.ind, J)
  } else {  # don't fix the number of individuals in each community
    n.ind <- round(rnorm(J, n.ind, 10))  
    n.ind[n.ind <= 0] <- n.ind  # set to n.ind if any generated number less than 1
  }
  full.data <- NULL
  
  for(j in 1:J) {
    if (verbose) message("#### generating " %+% j %+% "th cluster ####")
    cluster.data.j <- get.cluster.dat.Acont(id = j, n.ind = n.ind[j], truncBD = truncBD, 
                                            shift.val = shift.val, rndseed = eval(rndseed + j)) 
    full.data <- rbind(full.data, cluster.data.j)
  }  
  full.data$id <- as.integer(full.data$id)
  return(full.data)
}

`%+%` <- function(a, b) paste0(a, b)
J <- 1000
n.ind <- 100
rndseed <- 12345
truncBD <- 5
shift.val <- 1
comSample.wmT.cA.bY <- get.fullDat.Acont(J = J, n.ind = n.ind, rndseed = rndseed, 
                                         truncBD = truncBD, shift.val = shift.val)
N <- NROW(comSample.wmT.cA.bY)
mean(comSample.wmT.cA.bY$Y)  # 0.3478304
mean(comSample.wmT.cA.bY$Y.gstar)  # 0.4642384

define_f.gstar <- function(shift.val, truncBD, rndseed = NULL) {
  shift.const <- shift.val
  trunc.const <- truncBD
  f.gstar <- function(data, ...) {
    print("shift.const: " %+% shift.const)
    set.seed(rndseed)
    A.mu <- - 1.2 + 0.8 * data[, "E1"] + 0.21 * data[,"E2"] + 3 * mean(data[,"W2"]) - 
      0.7 * mean(data[,"W2"]) + 0.3 * mean(data[,"W3"])
    untrunc.A <- rnorm(n = nrow(data), mean = A.mu + shift.const, sd = 1)
    r.new.A <- exp(1.5 * shift.const * (untrunc.A - A.mu - shift.const / 4))
    trunc.A <- ifelse(r.new.A > trunc.const, untrunc.A - shift.const, untrunc.A)
    return(trunc.A)
  }
  return(f.gstar)
}
f.gstar <- define_f.gstar(shift.val = shift.val, truncBD = truncBD)

#*************************************** 
## Test 2.1 Different analysis methods
#*************************************** 
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)

## Test 2.1.1 Community-level analysis without a pooled individual-level regression on outcome
test_that("fit TMLE for cont, community-level A, with community-level analysis without a pooled Q", {
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4519760, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4498894, tolerance = 0.02)
  expect_equal(estimates["gcomp", ], 0.4567250, tolerance = 0.02)
})

## Test 2.1.2 Community-level analysis with a pooled individual-level regression on outcome
test_that("fit TMLE for cont, community-level A, with community-level analysis with a pooled Q", {
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = TRUE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4462443, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4510168, tolerance = 0.02)
  expect_equal(estimates["gcomp", ], 0.4494338, tolerance = 0.02)
})

# Test 2.1.3 Individual-level analysis with both individual-level outcome and treatment mechanisms
test_that("fit TMLE for cont, community-level A, with individual-level analysis", {
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "individual_level", communityID = "id",
                              rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4497872, tolerance = 0.01)
  expect_equal(estimates["iptw", ], 0.4519750, tolerance = 0.01)
  expect_equal(estimates["gcomp", ], 0.4525564, tolerance = 0.01)
})

# Test 2.1.4 Individual-level analysis that assumes no hierarchical structure
test_that("Treat hierarchical data as non-hierarchical and use usual TMLE", {
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "NoCommunity", rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4503520, tolerance = 0.01)
  expect_equal(estimates["iptw", ], 0.4526493, tolerance = 0.01)
  expect_equal(estimates["gcomp", ], 0.4531318, tolerance = 0.01)
  
  expect_message(
    tmleCom_res2 <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                                 WEnodes = c("E1", "E2", "W1", "W2", "W3"), 
                                 f_gstar1 = f.gstar, community.step = "individual_level", 
                                 communityID = NULL, rndseed = 12345),
    regexp = "Lack of 'communityID' forces the algorithm to automatically " %+%
      "pool data over all communities and treat it as non-hierarchical dataset"
  )
  estimates2 <- tmleCom_res2$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates, estimates2, tolerance = 0.01)
})

#*************************************** 
## Test 2.2 Different number of bins
#***************************************

## Test 2.2.1 Community-level analysis without a pooled Q + 10 bins
test_that("For cont, community-level A, community-level analysis without a pooled Q + 10 bins", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                  maxNperBin = N, nbins = 10)
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4504501, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4513750, tolerance = 0.02)
  expect_equal(estimates["gcomp", ], 0.4568511, tolerance = 0.02)
})

## Test 2.2.2 Community-level analysis without a pooled Q + 20 bins
test_that("For cont, community-level A, community-level analysis without a pooled Q + 20 bins", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                  maxNperBin = N, nbins = 50)
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4502626, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4689095, tolerance = 0.05)
  expect_equal(estimates["gcomp", ], 0.4563251, tolerance = 0.02)
})

#*************************************** 
## Test 2.3 Different binarization methods
#***************************************
## Test 2.3.1 Community-level analysis without a pooled Q + equal length
test_that("For cont, community-level A, community-level analysis without a pooled Q + equal length", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                  maxNperBin = N, bin.method = "equal.len")
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4503534, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4125413, tolerance = 0.02)
  expect_equal(estimates["gcomp", ], 0.4520120, tolerance = 0.02)
})

## Test 2.3.2 Community-level analysis without a pooled Q + combination of equal length & mass
test_that("For cont, community-level A, community-level analysis without a pooled Q + dhist", {
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", 
                  maxNperBin = N, bin.method = "dhist")
  tmleCom_res <-tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                              WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                              community.step = "community_level", communityID = "id", 
                              pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258
  expect_equal(estimates["tmle", ], 0.4521775, tolerance = 0.02)
  expect_equal(estimates["iptw", ], 0.4346268, tolerance = 0.02)
  expect_equal(estimates["gcomp", ], 0.4577096, tolerance = 0.02)
})

#*************************************** 
## Test 2.4 Different ML packages
#***************************************

## Test 2.4.1 Super Learner
test_that("fit TMLE for cont, community-level A with SL, using SL.glm, SL.bayesglm, SL.gam", {
  require("SuperLearner")
  tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = N,
                  SL.library = c("SL.glm", "SL.bayesglm", "SL.gam"), nbins = 5)
  tmleCom_res <- tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                               WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                               community.step = "community_level", communityID = "id", 
                               pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258 
  expect_equal(estimates["tmle", ], 0.4522367, tolerance = 0.02) 
  expect_equal(estimates["iptw", ], 0.4469491, tolerance = 0.02)  
  expect_equal(estimates["gcomp", ], 0.4554923, tolerance = 0.02) 
})

## Test 2.4.2 sl3
test_that("fit TMLE for cont, community-level A with sl3, using Lrnr_glm_fast", {
  require("sl3"); require("SuperLearner")
  tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "sl3_pipelines", maxNperBin = N, nbins = 5,
                  sl3_learner = list(glm_fast = sl3::make_learner(sl3::Lrnr_glm_fast)), 
                  sl3_metalearner = sl3::make_learner(sl3::Lrnr_optim, loss_function = sl3::loss_loglik_binomial,
                                                      learner_function = sl3::metalearner_logistic_binomial))
  tmleCom_res <- tmleCommunity(data = comSample.wmT.cA.bY, Ynode = "Y", Anodes = "A",
                               WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = f.gstar, 
                               community.step = "community_level", communityID = "id", 
                               pooled.Q = FALSE, rndseed = 12345)
  estimates <- tmleCom_res$EY_gstar1$estimates  # psi0 = 0.4571258 
  expect_equal(estimates["tmle", ], 0.4522367, tolerance = 0.02) 
  #expect_equal(estimates["iptw", ], 0.8559739, tolerance = 0.02)  # Need further investigation
  expect_equal(estimates["gcomp", ], 0.4554923, tolerance = 0.02) 
})
