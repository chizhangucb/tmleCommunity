context("Test tmleCommunity on individual-level exposure where outcome is rare (AKA CASE-CONTROL STUDY)")

`%+%` <- function(a, b) paste0(a, b)
Qform.corr <- "Y ~ W1 + W2*A + W3 + W4" # correct Q form
Qform.mis <- "Y ~ W3 + A" # incorrect Q form
gform.corr <- "A ~ W1 + W2 + W3 + W4"  # correct g form
gform.mis <- "A ~ W2"  # incorrect g form

# ---------------------------------------------------------------------------------
# TEST SET 1. Sample set with J = 1 (i.e., nCase/nControl = 1)
# ---------------------------------------------------------------------------------
data("indSample.iid.bA.bY.rareJ1_list", package = "tmleCommunity")
indSample.iid.bA.bY.rareJ1 <- indSample.iid.bA.bY.rareJ1_list$indSample.iid.bA.bY.rareJ1 
obs.wt.J1 <- indSample.iid.bA.bY.rareJ1_list$obs.wt.J1
N <- NROW(indSample.iid.bA.bY.rareJ1)
q0 <- indSample.iid.bA.bY.rareJ1_list$q0  # 0.013579
psi0.Y <- indSample.iid.bA.bY.rareJ1_list$psi0.Y  # 0.012662

test_that("Using correct observation weights and correctly specified Qform & gform, with J = 1", {
  tmleCom_Options()
  tmleCom_res <-
    tmleCommunity(data = indSample.iid.bA.bY.rareJ1, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                  obs.wts = obs.wt.J1)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.01220298, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.01277040, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.01243220, tolerance = 1e-6) 
})

test_that("Failing to specify obs weights and using correctly specified Qform & gform, with J = 1", {
  tmleCom_Options()
  tmleCom_res <-
    tmleCommunity(data = indSample.iid.bA.bY.rareJ1, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.2466575, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.2688000, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.2244907, tolerance = 1e-6) 
  
  tmleCom_res2 <-
    tmleCommunity(data = indSample.iid.bA.bY.rareJ1, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                  obs.wts = NULL)
  estimates2 <- tmleCom_res2$ATE$estimates  
  expect_equal(estimates, estimates2)
})

# ---------------------------------------------------------------------------------
# TEST SET 2. Sample set with J = 2 (i.e., nCase/nControl = 1/2)
# ---------------------------------------------------------------------------------
data("indSample.iid.bA.bY.rareJ2_list", package = "tmleCommunity")
indSample.iid.bA.bY.rareJ2 <- indSample.iid.bA.bY.rareJ2_list$indSample.iid.bA.bY.rareJ2 
obs.wt.J2 <- indSample.iid.bA.bY.rareJ2_list$obs.wt.J2
N <- NROW(indSample.iid.bA.bY.rareJ2)
psi0.Y <- indSample.iid.bA.bY.rareJ1_list$psi0.Y  # 0.012662

test_that("Using correct observation weights and correctly specified Qform & gform, with J = 2", {
  tmleCom_Options()
  tmleCom_res <-
    tmleCommunity(data = indSample.iid.bA.bY.rareJ2, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr,
                  obs.wts = obs.wt.J2)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.01228557, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.01273315, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.01234571, tolerance = 1e-6) 
})

test_that("Using correct obs weights & correctly specified gform & misspecified Qform, with J = 2", {
  tmleCom_Options()
  tmleCom_res <-
    tmleCommunity(data = indSample.iid.bA.bY.rareJ2, Ynode = "Y", Anodes = "A", 
                  WEnodes = c("W1", "W2", "W3", "W4"), f_gstar1 = 1, f_gstar2 = 0, 
                  Qform = Qform.mis, obs.wts = obs.wt.J2)
  estimates <- tmleCom_res$ATE$estimates  # psi0 = 0.012662 
  expect_equal(estimates["tmle", ], 0.012707583, tolerance = 1e-6)  
  expect_equal(estimates["iptw", ], 0.01273315, tolerance = 1e-6)  
  expect_equal(estimates["gcomp", ], 0.009430687, tolerance = 1e-6) 
})
