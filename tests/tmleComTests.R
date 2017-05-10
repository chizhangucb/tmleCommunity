load("Case I/EDA/data.Rda")
YsName <- c("premature_est", "zscore_bwgtgest")
AsName <- c("viol_r_tri1", "viol_r_tri2")
covEfix.Name <- c("pct_foreign_born", "pct_poverty", "pct_black", "pct_hisp", "pct_renters")
covEextra.Name <- c("umemp_a", "maxt_a", "pcp_a")
covEName <- c(covEfix.Name, covEextra.Name)
covWName <- c(c(paste0("mrace", 1:5), paste0("meduc_cat", 1:4), paste0("month", 1:11), 
                paste0("yr", 2007:2011)), "age_mom", "parity", "pay_gov_ins", "foreign_born")

create_f.gstar <- function(shiftPercent, truncBD = NULL, truncFUN = NULL, covNames) {
  shift.percent <- shiftPercent
  if (is.function(truncFUN)) { 
    truncFUN <- truncFUN
  } else if (all(is.number(truncBD))) {
    trunc.bound <- truncBD
  } else {
    stop("Either truncBD or truncFUN should be defined," %+%
           "where truncBD is a list of numbers and truncFUN is a function")
  }
  trtsName <- covNames
  f.gstar <- function(data, ...) {
    if (is.null(truncBD)) trunc.bound <- c(truncFUN(data[, trtsName[1]]), truncFUN(data[, trtsName[2]]))
    print("shift percentage: " %+% shift.percent)
    print("truncated bound for " %+% trtsName[1] %+% ": " %+% trunc.bound[1])
    print("truncated bound for " %+% trtsName[2] %+% ": " %+% trunc.bound[1])
    untrunc.A1 <- data[, trtsName[1]] * shift.percent
    untrunc.A2 <- data[, trtsName[2]] * shift.percent
    trunc.A1 <- base:::ifelse(trunc.bound[1] > untrunc.A1, trunc.bound[1], untrunc.A1)
    trunc.A2 <- base:::ifelse(trunc.bound[2] > untrunc.A2, trunc.bound[2], untrunc.A2)
    return(cbind(trunc.A1, trunc.A2))
  }
  return(f.gstar)
}
f.gstar <- create_f.gstar(shiftPercent = 0.8, truncFUN = min, covNames = AsName)

data <- subdata
gvars$verbose <- TRUE

######################################### 
## Test 1.1 speed.glm & glm 
######################################### 
tmleCom_Options(maxNperBin = nrow(data))
tmleCommunity_res.glm1 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = f.gstar, Qform = NULL, n_MCsims = 1)

tmleCommunity_res.glm1$EY_gstar1$estimates
#        estimate
# tmle   0.07298817
# h_iptw 0.07444963
# gcomp  0.07021494

tmleCommunity_res.glm2 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCommunity_res.glm2$EY_gstar1$estimates
#        estimate
# tmle    0.07016
# h_iptw  0.07016
# gcomp   0.07016
tmleCommunity_res.glm2$EY_gstar2$estimates
#        estimate
# tmle   0.07298817
# h_iptw 0.07444963
# gcomp  0.07021494
tmleCommunity_res.glm2$ATE$estimates
#             estimate
# tmle   -2.828174e-03
# h_iptw -4.289633e-03
# gcomp  -5.494195e-05


######################################### 
## Test 1.2 h2o
#########################################
tmleCom_Options(estimator = "h2o__ensemble", h2olearner =  c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), maxNperBin = nrow(data))
tmleCom_res.h2o.1 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                    Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCom_res.h2o.1$EY_gstar1$estimates
#          estimate
# tmle   0.07016000
# h_iptw 0.07016000
# gcomp  0.07070904
tmleCom_res.h2o.1$EY_gstar2$estimates
#           estimate
# tmle   0.071036598
# h_iptw 0.009401554
# gcomp  0.070604479
tmleCom_res.h2o.1$ATE$estimates
#             estimate
# tmle   -0.0008765975
# h_iptw  0.0607584455
# gcomp   0.0001045629

tmleCom_Options(estimator = "h2o__ensemble", maxNperBin = nrow(data),
                h2olearner =  as.vector(paste('h2o.glm.', c(0, 0.50, 1), sep = "")))
tmleCom_res.h2o.2 <- tmleCommunity(data = data, Ynode = YsName[1], Anodes = AsName, Wnodes = covWName, 
                                   Enodes = covEName, f_gstar1 = NULL, f_gstar2 = f.gstar, n_MCsims = 1)
tmleCom_res.h2o.2$EY_gstar1$estimates
#          estimate
# tmle   0.0701600
# h_iptw 0.07016000
# gcomp  0.07011594
tmleCom_res.h2o.2$EY_gstar2$estimates
#          estimate
# tmle   0.07096359
# h_iptw 0.08344168
# gcomp  0.07021925
tmleCom_res.h2o.2$ATE$estimates
#             estimate
# tmle   -0.0008035862
# h_iptw -0.0132816815
# gcomp  -0.0001033166

