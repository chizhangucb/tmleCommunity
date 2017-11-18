context("Test GenericModel")

library(data.table)
`%+%` <- function(a, b) paste0(a, b)
data("indSample.iid.cA.cY_list", package = "tmleCommunity")
indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
N <- nrow(indSample.iid.cA.cY)
tmleCom_Options(maxNperBin = N, nbins = 5, bin.method = "equal.mass")
nodes <- list(Ynode = "Y", Anodes = "A", WEnodes = c("W1", "W2", "W3", "W4"))
h.g0.sVars <- tmleCommunity:::define_regform(A ~ W1 + W2 + W3 + W4)
subsets_expr <- lapply(h.g0.sVars$outvars, function(var) {var})
OData.g0 <- DatKeepClass$new(Odata = indSample.iid.cA.cY, nodes = nodes)
regclass.g0 <- RegressionClass$new(outvar = h.g0.sVars$outvars,
                                   predvars = h.g0.sVars$predvars,
                                   subset_vars = subsets_expr,
                                   outvar.class = OData.g0$type.sVar[h.g0.sVars$outvars])
genericmodels.g0 <- GenericModel$new(reg = regclass.g0, DatKeepClass.g0 = OData.g0)
genericmodels.g0$fit(data = OData.g0)
genericmodels.g0.A1 <- genericmodels.g0$getPsAsW.models()[[1]]

test_that("Different entrances to fit regression bin(s) and get the same results", {
  # Fit GenericModel directly through the parent RegressionClass
  genericmodels.g0.A1.B2 <- genericmodels.g0.A1$getPsAsW.models()[[2]]
  # genericmodels.g0.A1.B2$getfit$coef
  
  # Fit ContinModel through the clone of a parent RegressionClass 
  k_i <- 1  # The first exposure node
  reg_i <- regclass.g0$clone()
  reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
  genericmodels.g0.A1.cont <- ContinModel$new(reg = reg_i, DatKeepClass.g0 = OData.g0)
  genericmodels.g0.A1.cont$fit(data = OData.g0, savespace = FALSE)
  genericmodels.g0.A1.cont.B2 <- genericmodels.g0.A1.cont$getPsAsW.models()[[2]]
  # genericmodels.g0.A1.cont.B2$getfit$coef
  
  # Fit GenericModel directly from genericmodels.g0.A1.alt (The first exposure node)
  subset_eval.A1 <- tmleCommunity:::def_regs_subset(genericmodels.g0.A1.cont)
  genericmodels.g0.A1.cont.Gen <- GenericModel$new(reg = subset_eval.A1)
  genericmodels.g0.A1.cont.Gen.B2 <- genericmodels.g0.A1.cont.Gen$getPsAsW.models()[[2]]
  genericmodels.g0.A1.cont.Gen.B2$fit(data = OData.g0)
  genericmodels.g0.A1.cont.Gen.B2$getfit$coef 
  
  expect_equal(genericmodels.g0.A1.B2$getfit$coef, genericmodels.g0.A1.cont.B2$getfit$coef)
  expect_equal(genericmodels.g0.A1.B2$getfit$coef, genericmodels.g0.A1.cont.Gen.B2$getfit$coef)
})

test_that("The first bin contains all zeros and so produce tiny coefficients; " %+%  
            "The last bin should contain nothing and so produce NA coefficients", {
              lastB.Ind <- eval(regclass.g0$nbins + 2)          
              genericmodels.g0.A1.B1 <- genericmodels.g0.A1$getPsAsW.models()[[1]]
              genericmodels.g0.A1.Bend <- genericmodels.g0.A1$getPsAsW.models()[[lastB.Ind]]          
              expect_true(all(genericmodels.g0.A1.B1$getfit$coef[-1] < 10e-8))
              expect_true(all(is.na(genericmodels.g0.A1.Bend$getfit$coef)))          
            })

test_that("Binarized bins contain the (approximately) same number of obs if 'equal.mass'", {
  OData.g0$binirize.sVar(name.sVar = genericmodels.g0.A1$outvar, intervals = genericmodels.g0.A1$intrvls, 
                         nbins = genericmodels.g0.A1$nbins, bin.nms = genericmodels.g0.A1$bin_nms)
  ord.sVar <- tmleCommunity:::make.ordinal(x = OData.g0$get.sVar(genericmodels.g0.A1$outvar), 
                                           intervals = genericmodels.g0.A1$intrvls)
  expect_true(all(table(ord.sVar) == N / regclass.g0$nbins))
  expect_equal(sum(ord.sVar == 1), 0)  # NO obs falls beyond the minimum bound
  expect_equal(sum(ord.sVar == eval(regclass.g0$nbins + 2)), 0)  # NO obs falls beyond the max bound
})

test_that("Pool bin indicators across all bins and fit one pooled regression", {
  k_i <- 1  # The first exposure node
  reg_i <- regclass.g0$clone()
  reg_i$ChangeManyToOneRegresssion(k_i, regclass.g0)
  datX_mat <- OData.g0$get.dat.sVar(rowsubset = TRUE, covars = (genericmodels.g0.A1$bin_nms))
  pooled_bin_name <- OData.g0$pooled.bin.nm.sVar(reg_i$outvar)  # "A_allB.j"
  expect_warning(  # Expect no warning when pooling across all bins
    BinsDat_long <- 
      tmleCommunity:::binirized.to.DTlong(BinsDat_wide = datX_mat, binID_seq = (1L:eval(reg_i$nbins + 2)), 
                                          ID = as.integer(1:nrow(datX_mat)), 
                                          bin_names = genericmodels.g0.A1$bin_names, 
                                          pooled_bin_name = pooled_bin_name, name.sVar = reg_i$outvar),
    regexp = NA)
  sVar_melt_DT <- 
    tmleCommunity:::join.Xmat(X_mat = OData.g0$get.dat.sVar(TRUE, covars = reg_i$predvars), 
                              sVar_melt_DT = BinsDat_long, ID = as.integer(1:nrow(datX_mat)))
  # select bin_ID + predictors, add intercept column
  X_mat <- sVar_melt_DT[,c("bin_ID", reg_i$predvars), with=FALSE][, c("Intercept") := 1]
  # re-order columns by reference (no copy)
  setcolorder(X_mat, c("Intercept", "bin_ID", reg_i$predvars)) 
  Y_vals <- sVar_melt_DT[, pooled_bin_name, with = FALSE][[1]]
  expect_equal(names(X_mat), c("Intercept", "bin_ID", nodes$WEnodes))
  expect_true(all(X_mat$bin_ID %in% (1 : eval(reg_i$nbins + 2))))
  # Only the bin where the obs falls into receives 1, others receive 0
  expect_equal(sort(unique(Y_vals)), c(0, 1))  
})

test_that("Test if predictAeqa() provides the same results as the cumulative prob * bin weights", {
  intrvls.width <- diff(genericmodels.g0.A1$intrvls)
  bws <- OData.g0$get.sVar.bw(name.sVar = genericmodels.g0.A1$outvar, 
                              intervals = genericmodels.g0.A1$intrvls)
  bin_weights <- (1 / bws)  
  genericmodels.g0.A1.B2 <- genericmodels.g0.A1$getPsAsW.models()[[2]]
  genericmodels.g0.A1.B2$newdata(newdata = OData.g0, getoutvar = TRUE)
  probA1 <- tmleCommunity:::predict_single_reg(self = genericmodels.g0.A1.B2)
  probA1 <- probA1[genericmodels.g0.A1.B2$subset_idx]
  indA <- OData.g0$get.outvar(genericmodels.g0.A1.B2$subset_idx, var = genericmodels.g0.A1.B2$outvar)
  probAeqa <- rep.int(1L, OData.g0$nobs) 
  probAeqa[genericmodels.g0.A1.B2$subset_idx] <- probA1^(indA) * (1 - probA1)^(1L - indA)
  genericmodels.g0.A1$predictAeqa(newdata = OData.g0, wipeProb = FALSE)
  expect_equal(mean(probAeqa), mean(genericmodels.g0.A1.B2$getprobAeqa))
  
  cumprodAeqa <- 1
  for (i in 1:length(genericmodels.g0.A1$getPsAsW.models())) {
    cumprodAeqa <- cumprodAeqa * genericmodels.g0.A1$getPsAsW.models()[[i]]$getprobAeqa
  }
  expect_equal(mean(cumprodAeqa * bin_weights), mean(genericmodels.g0.A1$getcumprodAeqa()))
})
