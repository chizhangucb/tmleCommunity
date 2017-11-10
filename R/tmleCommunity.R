# @title tmleCommunity-package
# @docType package
# @author Chi Zhang, Oleg Sofrygin, Mark J. van der Laan

#' @import R6
#' @importFrom grDevices nclass.FD nclass.Sturges nclass.scott
#' @importFrom graphics axis barplot hist par text
#' @importFrom methods is
#' @importFrom stats approx quasibinomial binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm runif 
#' @importFrom stats terms var aggregate weighted.mean as.formula pnorm pt qt reformulate model.matrix gaussian 
#' @importFrom utils data head str
#' @importFrom Hmisc wtd.var
#' @importFrom assertthat assert_that is.flag is.string is.count

# ------------------------------------ define_regform ----------------------------------------- 
# Purpose: Parse the formulae for actual covariate names in A & W
# --------------------------------------------------------------------------------------------- 
define_regform <- function(regform, Anodes.lst = NULL, Wnodes.lst = NULL) {
  if (length(regform)==0L) {
    return(list(outvars =  as.vector(unlist(Anodes.lst)), predvars = as.vector(unlist(Wnodes.lst))))
  } else {
    # Getting predictors (W names):
    regformterms <- terms(regform)
    W.names <- attributes(regformterms)$term.labels 
    W.names.alt <- colnames(attributes(regformterms)$factors)
    assert_that(all(W.names == W.names.alt))
    
    # Getting outcomes (A names):
    out.var <- rownames(attributes(regformterms)$factors)[1] # character string
    out.vars.form <- as.formula(". ~ " %+% out.var)
    out.vars.terms <- terms(out.vars.form)
    A.names <- attributes(out.vars.terms)$term.labels
    
    # (Interface for specifying regressions for g0_form & gstar_form)
    get_vars_fromlist <- function(varname, Var.lst) {
      if (varname %in% names(Var.lst)) { as.vector(Var.lst[[varname]]) } else { varname }
    }
    
    outvars <- unlist(lapply(A.names, get_vars_fromlist, Anodes.lst))
    predvars <- unlist(lapply(W.names, get_vars_fromlist, Wnodes.lst))
    return(list(outvars = outvars, predvars = predvars))
  }
}


#--------------------------------------- tmle.update ---------------------------------------
# Purpose: updating estiamted Qbar by clever covariate
#-------------------------------------------------------------------------------------------
tmle.update <- function(TMLE.targetStep, Y, obs.wts, off, h_wts, subset, family) {
  ctrl <- glm.control(trace = FALSE, maxit = 1000)
  if (family %in% c("binomial", "quasibinomial")) { glmFamily <- quasibinomial() } else { glmFamily <- gaussian() }
  if (TMLE.targetStep == "tmle.intercept") {
    # ************************************************
    # estimate the TMLE update via weighted univariate ML (espsilon is intercept)
    # ************************************************
    # SuppressGivenWarnings(model.Q.star <- glm(Y ~ offset(off), data = data.frame(Y = Y, off = off), weights = h_wts,
    #                                           subset = subset, family = glmFamily, control = ctrl), GetWarningsToSuppress(TRUE))
    SuppressGivenWarnings(
      model.Q.star <- stats::glm.fit(x = rep(1, length(Y))[subset], y = Y[subset], offset = off[subset], 
                                     family = glmFamily, weights = h_wts * obs.wts, control = ctrl), GetWarningsToSuppress(TRUE)
    )
    QY.star <- Y
    if (!is.na(coef(model.Q.star))) QY.star <- plogis(off + coef(model.Q.star))
  } else if (TMLE.targetStep == "tmle.covariate") {
    # ************************************************
    # estimate the TMLE update via univariate ML (epsilon is coefficient for h^*/h) - ONLY FOR NON-DETERMINISTIC SUBSET   
    # ************************************************
    # SuppressGivenWarnings(model.Q.star <- glm(Y ~ -1 + h_wts + offset(off), data = data.frame(Y = Y, off = off, h_wts = h_wts),
    #                                           subset = subset, family = glmFamily, control = ctrl), GetWarningsToSuppress(TRUE))
    SuppressGivenWarnings(
      model.Q.star <- stats::glm.fit(x = as.matrix(data.frame(h_wts = h_wts))[subset, ], y = Y[subset], offset = off[subset], 
                                     weights = obs.wts, family = glmFamily, control = ctrl), GetWarningsToSuppress(TRUE)
    )
    QY.star <- Y
    if (!is.na(coef(model.Q.star))) QY.star <- plogis(off + coef(model.Q.star) * h_wts)
  }
  return(list(model.Q.star = model.Q.star, QY.star = QY.star))
}


#------------------------------------- calcParameters --------------------------------------
# Purpose: create output object with param ests of EY_gstar, vars and CIs for given gstar 
# (or ATE if two tmle obj are passed) boot.var, n.boot
#-------------------------------------------------------------------------------------------
calcParameters <- function(inputYs, alpha = 0.05, est_params_list, tmle_g_out, tmle_g2_out = NULL) {
  OData.ObsP0 <- tmle_g_out$OData.ObsP0
  # If "perCommunity", still use the number of independent communities in variance calculation (aggregate IC by community first)
  if (est_params_list$community.step %in% c("community_level", "individual_level", "perCommunity")) {
    est_params_list$communityID <- OData.ObsP0$get.sVar(name.sVar = est_params_list$communityID)
    nobs <- length(unique(est_params_list$communityID))
  } else {
    nobs <- OData.ObsP0$nobs
  }
  df <- ifelse(nobs <= 40, (nobs - 2), NA)  # Use the Student's T distribution in place of the Std Normal if nobs < 41
  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat
  obs.wts <- tmle_g_out$obs.wts
  maptoYstar <- inputYs$maptoYstar
  ab <- inputYs$ab

  if (!is.null(tmle_g2_out)) {
    ests_mat <- tmle_g_out$ests_mat - tmle_g2_out$ests_mat
    fWi_mat <- tmle_g_out$fWi_mat - tmle_g2_out$fWi_mat
    wts_mat <- tmle_g_out$wts_mat - tmle_g2_out$wts_mat
  }
  
  # *****************************************************
  # get the iid IC-based asymptotic variance estimates:
  # *****************************************************
  var_mat.res <- get_est_sigmas(estnames = c("tmle", "iptw", "gcomp"), obsYvals = OData.ObsP0$noNA.Ynodevals, est_params_list= est_params_list, 
                                obs.wts = obs.wts, ests_mat = ests_mat, QY_mat = QY_mat, wts_mat = wts_mat, fWi_mat = fWi_mat)
  as.var_mat <- var_mat.res$as.var_mat
  if (maptoYstar) {
    as.var_mat <- as.var_mat * (diff(ab) ^ 2)
    if (is.null(tmle_g2_out)) {
      ests_mat <- ests_mat * diff(ab) + ab[1]
    } else { 
      ests_mat <- ests_mat * diff(ab) 
    }
  }
  
  # Test statistic & p-value
  teststat_mat <-  matrix(0L, nrow = 3, ncol = 1)
  teststat_mat[, 1] <- ests_mat[, 1] / sqrt(as.var_mat / nobs)
  if (is.na(df)) {
    pval <- 2 * pnorm(abs(teststat_mat), lower.tail = F) 
  } else {
    pval<- 2 * pt(abs(teststat_mat), df = df, lower.tail = F) 
  }

  get_CI <- function(xrow, n, df = NA) {
    f_est_CI <- function(n, psi, sigma2_N, df) { # get CI
      if (is.na(df)) {
        cutoff <- qnorm(1 - alpha/2) # z_alpha 
      } else {
        cutoff <- qt(1 - alpha/2, df = df) # t_alpha 
      }
      CI_est <- c(psi - cutoff * sqrt(sigma2_N / n), psi + cutoff * sqrt(sigma2_N / n))
      return(CI_est)
    }
    psi <- xrow["estimate"];
    sigma2_N <- xrow["Var"];
    return(f_est_CI(n = n, psi = psi, sigma2_N = sigma2_N, df = df))
  }
  
  CIs_mat <- t(apply(cbind(ests_mat, as.var_mat), 1, get_CI, n = nobs, df = df))
  colnames(CIs_mat) <- c("LBCI_" %+% as.character(alpha/2), "UBCI_" %+% as.character(1-alpha/2))
  
  # ----------------------------------------------------
  # RENAME ESTIMATORS FOR THE FINAL OUTPUT:
  # ----------------------------------------------------
  rownames(ests_mat) <- c("tmle", "iptw", "gcomp")
  rownames(as.var_mat) <- c("tmle", "iptw", "gcomp")
  rownames(CIs_mat) <- c("tmle", "iptw", "gcomp")
  rownames(teststat_mat) <- rownames(pval) <- c("tmle", "iptw", "gcomp")
  colnames(teststat_mat) <- "teststat"; colnames(pval) <- "p_value"
  
  EY_g.star <- list(estimates = ests_mat, 
                    vars = (as.var_mat / nobs), 
                    CIs = CIs_mat, 
                    tstat = teststat_mat,
                    pval = pval,
                    IC = var_mat.res$IC, 
                    h.g0_GenericModel = tmle_g_out$h.g0_GenericModel, 
                    h.gstar_GenericModel = tmle_g_out$h.gstar_GenericModel)
  
  if (!is.null(tmle_g2_out)) {
    EY_g.star[["h.g0_GenericModel"]] <- NULL
    EY_g.star[["h.gstar_GenericModel"]] <- NULL
  }
  return(EY_g.star)
}


#------------------------------------------ get_est_sigmas ---------------------------------------
# Purpose: get the iid IC-based asymptotic variance estimates
# Formula source: http://biostats.bepress.com/cgi/viewcontent.cgi?article=1351&context=ucbbiostat (Page 18)
# OR http://biostats.bepress.com/cgi/viewcontent.cgi?article=1292&context=ucbbiostat (Page 5)
#-------------------------------------------------------------------------------------------------
get_est_sigmas <- function(estnames, obsYvals, est_params_list, obs.wts, ests_mat, QY_mat, wts_mat, fWi_mat) {
  community.step <- est_params_list$community.step
  community.wts <- est_params_list$community.wts
  communityID <- est_params_list$communityID
  
  fWi <- fWi_mat[, "fWi_Qinit"]
  QY.init <- QY_mat[, "QY.init"] 
  h_wts <- wts_mat[, "h_wts"]
  
  # TMLE inference based on the iid IC: (** Use QY.init not QY.star)
  iidIC_tmle <- h_wts * (obsYvals - QY.init) + fWi - ests_mat[rownames(ests_mat) %in% "TMLE",]
  # MLE inference based on the iid IC: (** Use QY.init not QY.star) *** NOT ACCURATE ***
  iidIC_mle <- h_wts * (obsYvals - QY.init) + fWi - ests_mat[rownames(ests_mat) %in% "MLE",]
  # IPTW h (based on the mixture density clever covariate (h)):
  iidIC_iptw <- h_wts * (obsYvals) - ests_mat[rownames(ests_mat) %in% "IPTW",]
  
  iidIC <- data.table::data.table(iidIC_tmle, iidIC_mle, iidIC_iptw, obs.wts); sorted.communityID <- communityID
  # if we believe our working model (i.e. if estimating under the submodel) or run TMLE for each community
  if (community.step == "individual_level" || community.step == "perCommunity") { 
    # if (!is.null(communityID)) {"iid IC cannnot be aggregated to the cluster-level since lack of 'communityID' so treated as non-hierarchical"}
    iidIC <- data.table::data.table(iidIC_tmle, iidIC_mle, iidIC_iptw, obs.wts, communityID)
    iidIC <- iidIC[, lapply(.SD, weighted.mean, w = obs.wts), by = communityID]; sorted.communityID <- iidIC[["communityID"]]
    obs.wts <- community.wts[match(sorted.communityID, community.wts[, "id"]), "weights"]
    # iidIC_tmle <- aggregate(x = iidIC_tmle, by=list(newid = communityID), mean), similarly to iidIC_mle, iidIC_iptw
    # obs.wts <- community.wts[match(iidIC_tmle[, 1], community.wts[, "id"]), "weights"]  # Alternative way
  }
  
  iidIC <- iidIC[, !(colnames(iidIC) %in% c("obs.wts", "communityID")), with = FALSE]
  var_iid.est <- iidIC[, lapply(.SD, Hmisc::wtd.var, weights = obs.wts, normwt = T)]
  # var_iid.tmle <- Hmisc::wtd.var(iidIC_tmle, weights = obs.wts, normwt = T), similarly to var_iid.mle, var_iid.iptw
  # var_iid.tmle <- mean((iidIC_tmle)^2)  # the same as sum(iidIC_tmle^2) / length(iidIC_tmle) by assume mean(iidIC_tmle) = 0
  as.var_mat <- matrix(0, nrow = 3, ncol = 1)
  as.var_mat[, 1] <- c(var_iid.est[["iidIC_tmle"]], var_iid.est[["iidIC_iptw"]], var_iid.est[["iidIC_mle"]])
  rownames(as.var_mat) <- estnames; colnames(as.var_mat) <- "Var"
  iidIC <- as.data.frame(iidIC); rownames(iidIC) <- sorted.communityID; colnames(iidIC) <- c("IC.tmle", "IC.iptw", "IC.gcomp")
  return(list(as.var_mat = as.var_mat, IC = iidIC))
}


#--------------------------------------- CalcAllEstimators ---------------------------------
# Purpose: Estimate h_bar under g0 and g* given observed data and vector of c^Y's data is an DatKeepClass object
#-------------------------------------------------------------------------------------------
CalcAllEstimators <- function(OData.ObsP0, est_params_list) {
  data <- est_params_list$data
  communityID <- est_params_list$communityID
  community.step <- est_params_list$community.step
  community.wts <- est_params_list$community.wts
  pooled.Q <- est_params_list$pooled.Q
  f.gstar <- est_params_list$f.gstar
  TMLE.targetStep <- est_params_list$TMLE.targetStep
  model.Q.init <- est_params_list$model.Q.init
  nodes <- OData.ObsP0$nodes
  Y <- OData.ObsP0$noNA.Ynodevals # actual observed & transformed Y's
  determ.Q <- OData.ObsP0$det.Y
  
  # getting all node vals, including deterministic Y and predictions P(Y=1) for non-DET Y under orignal A
  QY.init <- OData.ObsP0$noNA.Ynodevals 
  QY.init[!OData.ObsP0$det.Y] <- model.Q.init$predict(newdata = OData.ObsP0)$getprobA1[!OData.ObsP0$det.Y]
  QY.init <- bound(QY.init, est_params_list$Qbounds)
  off <- qlogis(QY.init)  # offset
  
  if (community.step == "community_level" && pooled.Q) { # if we do NOT believe our working model (i.e. estimate under the lareg model)
    # if (!is.null(communityID)) {}  # Don't need it here since tmleCommunity takes care of it
    # Since community-level TMLE with pooled individual-level Q requires 'communityID' to aggregate data to the cluster-level in the estimation  
    # of trt mechanism. Lack of 'communityID' pool data over all communities & treat it as non-hierarchical data when fitting clever covariates
    
    # aggregate initial outcome predictions to the cluster-level & Recalculate offset based on aggregated initiate predictions
    QY.init <- aggregate(x = QY.init, by=list(id = data[, communityID]), mean)[, 2]
    off <- qlogis(QY.init)  
    # aggregate original dataset to the cluster-level & redefine OData.ObsP0
    Y <- aggregate(x = Y, by=list(id = data[, communityID]), mean)[, 2] 
    determ.Q <- rep_len(FALSE, length(Y))  # For aggregated data, YnodeDet is currently unavailable, treat all Y^c as nondeterministic
    f_gstar.agg <- FALSE
    # Also need to aggregate f_gstar if it is a vector of length NROW(data) (or a matrix/ data.fram with nrow = NROW(data)) 
    if (is.VecMatDf(f.gstar) && NROW(as.data.frame(f.gstar)) == NROW(data)) { f_gstar.agg <- TRUE; data <- cbind(data, f.gstar) }
    data <- aggregate(x = data, by=list(newid = data[, communityID]), mean) 
    colname.allNA <- colnames(data)[colSums(is.na(data)) == NROW(data)]  # columns with all NAs after aggregation, due to non-numeric values
    if (length(colname.allNA) != 0) {
      data <- data[, - which(colnames(data) %in% colname.allNA)] # Remove columns that contain only NAs inside
      names(data)[which(colnames(data) == "newid")] <- communityID  # change 'newid' back to the communityID name
      warning(paste(colname.allNA, collapse = ', ') %+% " is(are) removed from the aggregated data due to all NAs in the column(s).")
      warning("Suggestion: convert the non-numeric values to numeric, e.g., create dummy variables for each category/ string/ factor.")
    }
    est_params_list$data <- data
    est_params_list$obs.wts <- obs.wts <- community.wts[match(data[, communityID], community.wts[, "id"]), "weights"]  # ensure matching
    if (f_gstar.agg) { 
      if (is.vector(f.gstar)) { est_params_list$f.gstar <- data$f.gstar } else {est_params_list$f.gstar <- data[, names(f.gstar)]} 
    } 
    OData.ObsP0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
    OData.ObsP0$addYnode(YnodeVals = data[, est_params_list$nodes$Ynode])  # Already bounded Y into Ystar in the beginning step               
    OData.ObsP0$addObsWeights(obs.wts = obs.wts)
  } 
  
  #************************************************
  # Fitting h_gstar / h_gN clever covariate:
  #************************************************
  fit.hbars_t <- system.time(fit.hbars.out <- fit.hbars(OData.ObsP0 = OData.ObsP0, est_params_list = est_params_list)) # fit the clever covariate
  OData.gstar <- fit.hbars.out$OData.gstar
  model.h.fit <- fit.hbars.out$model.h.fit
  h_wts <- fit.hbars.out$h_gstar_h_gN
  obs.wts <- est_params_list$obs.wts <- OData.gstar$get.obsweights(TRUE)  # Double check here
  
  #************************************************
  # IPTW_h estimator:
  #************************************************  
  IPTW <- Y
  IPTW[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]
  # IPTW_unwt <- mean(IPTW)
  if (community.step == "individual_level") {
    # if (!is.null(communityID)) {"IPTW cannnot be aggregated to the cluster-level since lack of 'communityID' so treated as non-hierarchical"}
    # IPTW <- aggregate(x = IPTW, by=list(newid = data[, communityID]), mean)
    newid <- data[, communityID]; IPTW <- cbind(IPTW, obs.wts)
    wts.mean <- function(d) { weighted.mean(x = d[, "IPTW"], w = d[, "obs.wts"]) }
    IPTW <- Hmisc::summarize(X = IPTW, by = newid, FUN = wts.mean)
    IPTW <- weighted.mean(IPTW[, 2], w = community.wts[match(IPTW[, "newid"], community.wts[, "id"]), "weights"])
  } else {
    IPTW <- weighted.mean(IPTW, w = obs.wts)
  }
  
  #************************************************
  # TMLE estimators
  #************************************************  
  tmle.update.out <- tmle.update(TMLE.targetStep = TMLE.targetStep, Y = Y, off = off, obs.wts = obs.wts,  
                                 h_wts = h_wts, subset = !determ.Q, family = "quasibinomial")
  model.Q.star <- tmle.update.out$model.Q.star
  QY.star <- tmle.update.out$QY.star
  
  #************************************************
  # Run Monte-Carlo (MC) evaluation for all plug-in estimators (TMLE & Gcomp), under stochastic intervention g^*:
  #************************************************
  MC_fit_params <- append(est_params_list, list(model.Q.star = model.Q.star))
  syst1 <- system.time(MCS_res <- CalcMonteCarloEsts(OData.ObsP0 = OData.ObsP0, 
                                                     OData.gstar = OData.gstar, 
                                                     MC_fit_params = MC_fit_params, 
                                                     model.h.fit = model.h.fit))

  ests <- c(TMLE = MCS_res[["TMLE"]], IPTW = IPTW, MLE = MCS_res[["MLE"]])
  ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
  ests_mat[, 1] <- ests
  rownames(ests_mat) <- names(ests); colnames(ests_mat) <- "estimate"
  
  wts_mat <- matrix(0L, nrow = OData.ObsP0$nobs, ncol = 1)
  colnames(wts_mat) <- c("h_wts")
  wts_mat[, "h_wts"] <- h_wts
  
  fWi_mat <- matrix(0L, nrow = OData.ObsP0$nobs, ncol = 1)
  colnames(fWi_mat) <- c("fWi_Qinit")
  fWi_mat[,"fWi_Qinit"] <- MCS_res[agrep("fWi_init_", names(MCS_res), max.distance=list(insertions=0, substitutions=0))]
  
  QY_mat <- matrix(0L, nrow = OData.ObsP0$nobs, ncol = 2)
  colnames(QY_mat) <- c("QY.init", "QY.star")
  QY_mat[,] <- cbind(QY.init, QY.star)
  
  if (gvars$verbose)  {
    print("time spent fitting new fit.hbars.out:"); print(fit.hbars_t)
    if (TMLE.targetStep == "tmle.covariate") {
      parsubmodel_fits <- rbind(coef(model.Q.star))
      rownames(parsubmodel_fits) <- c("epsilon (clever covariate coefficient)")
    } else if (TMLE.targetStep == "tmle.intercept") {
      parsubmodel_fits <- rbind(coef(model.Q.star))
      rownames(parsubmodel_fits) <- c("alpha (intercept)")
    }
    print("new parsubmodel_fits: "); print(parsubmodel_fits)
    print("time to run Monte Carlo target param evaluation: "); print(syst1);
    print(c(fWi_init = mean(fWi_mat[,"fWi_Qinit"] - ests["TMLE"])));
    print("new MC.ests mat: "); print(ests_mat)
  }
  
  return(list(ests_mat = ests_mat,
              wts_mat = wts_mat,
              fWi_mat = fWi_mat,
              QY_mat = QY_mat,
              obs.wts = obs.wts, 
              h.g0_GenericModel = model.h.fit$genericmodels.g0,
              h.gstar_GenericModel = model.h.fit$genericmodels.gstar,
              OData.ObsP0 = OData.ObsP0))
}


#------------------------------------
#' Estimate Marginal Treatment Effects For Arbitrary (Stochastic) Interventions in Hierarchical Data
#'
#' Estimate the marginal treatment effect among independent communities (or i.i.d units if no hierarchical structure) using \strong{TMLE} 
#' (targeted maximum likelihood estimation). It supports two different TMLEs that are based on community-level and individual-level 
#' analysis, respectively. The individual-level TMLE cooperates with additional working assumptions and has potential efficiency gain. It also  
#' provide corresponding \strong{IPTW} (the inverse-probability-of-treatment or Horvitz-Thompson) and \strong{GCOMP} (parametric G-computation).
#'
#' @param data Observed data, \code{data.frame} with named columns, containing \code{WEnodes}, \code{Anode}, \code{Ynode} and possibly 
#'  \code{communityID}, \code{YnodeDet}. See "Details".
#' @param Ynode Column names or indices in \code{data} of outcome variable. Outcome can be either binary or continuous (could be beyond 0 and 1). 
#'  If Ynode undefined, the left-side of the regression formula in argument \code{Qform} will be treated as \code{Ynode}.
#' @param Anodes Column names or indices in \code{data} of exposure (treatment) variables
#' @param WEnodes Column names or indices in \code{data} of individual-level (and possibly community-level) baseline covariates. Factors are not allowed.
#' @param YnodeDet Optional column name or index in \code{data} of indicators of deterministic values of outcome \code{Ynode}, coded as (TRUE / FALSE)  
#'  or (1 / 0). If TRUE or 1, value of \code{Ynode} is given deterministically / constant. 
#' @param obs.wts Optional choice to provide/ construct a vector of individual-level observation (sampling) weights (of length \code{nrow(data)}). 
#'  Currently supports a non-negative numeric vector, "equal.within.pop" (Default) and equal.within.community. If "equal.within.pop", weigh individuals 
#'  in the entire dataset equally (weigh to be all 1); If "equal.within.community", weigh individuals within the same community equally 
#'  (i.e., 1 / (number of individuals in each community)).
#' @param community.step Methods to deal with hierarchical data, one of "NoCommunity" (Default), "community_level", "individual_level" and "PerCommunity".  
#'  If "NoCommunity", claim that no hirerachical structure in data; If "community_level", use community-level TMLE; If "individual_level", use  
#'  individual-level TMLE cooperating with the assumption of no covariate interference. If "perCommunity", use stratified TMLE. If \code{communityID} =  
#'  \code{NULL}, then automatically pool over all communities (i.e., treated it as "NoCommunity"). See "Details".
#' @param communityID Optional column name or index in \code{data} of community identifier variable. If known, it can support the two options within 
#'  \code{community.step}: community-level or individual-level TMLE (See details for \code{community.step}).
#' @param community.wts Optional choice to provide/ construct a matrix of community-level observation weights (where dimension = \eqn{J\times}2, where 
#'  J = the number of communities). The first column contains the communities' names (ie., \code{data[, communityID]}) and the second column contains the   
#'  corresponding non-negative weights. Currently only support a numeric matrix with 2 columns, "size.community" (Default) and "equal.community". 
#'  If setting community.wts = "size.community", treat the number of individuals within each community as its weight, respectively. 
#'  If community.wts = "equal.community", assumed weights to be all 1; 
#' @param pooled.Q Logical for incorporating hierarchical data to estimate the outcome mechanism. If \code{TRUE}, use a pooled individual-level
#'  regression for initial estimation of the mean outcome (i.e., outcome mechanism). Default to be \code{FASLE}. See "Details". 
#' @param f_g0 Optional function used to specify model knowledge about value of Anodes. It estimates P(A | W, E) under \code{g0} by 
#'  sampling a large vector/ data frame of Anode (of length \code{nrow(data)*n_MCsims} or number of rows if a data frame) from \code{f_g0} function.
#' @param f_gstar1 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  If a matrix/ data frame, its number of rows must be either \code{nrow(data)} or 1 (constant exposure assigned to all observations), and its number of 
#'  columns must be \code{length(Anodes)}. Note that the column names should match with the names in \code{Anodes}. If a vector, it must be of length 
#'  \code{nrow(data)} or 1. If a function, it must return a vector or a data frame of counterfactual exposures sampled based on Anodes, WEnodes
#'  (and possibly communityID) passed as a named argument "data". Thus, the function must include "data" as one of its argument names. 
#'  The interventions defined by f_gstar1 can be static, dynamic or stochastic. See Exmaples below.
#' @param f_gstar2 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  It has the same components and requirements as f_gstar1 has.
#' @param Qform Character vector of regression formula for Ynode. If not specified (i.e., \code{NULL}), the outcome variable is regressed on all 
#'  covariates included in Anodes and WEnodes (i.e., \code{Ynode ~ Anodes + WEnodes}). See "Details".
#' @param Qbounds Vector of upper and lower bounds on Y and predicted value for initial Q. Default to the range of Y, widened by 10\% of the min 
#'  and max values. See "Details".
#' @param alpha Used to keep predicted values for initial Q bounded away from (0,1) for logistic fluctuation 
#'  (set \code{Qbounds} to (1 - \code{alpha}), \code{alpha}).
#' @param fluctuation Default to "logistic", it could also be "linear" (for targeting step).
#' @param hform.g0 Character vector of regression formula for estimating the conditional density of P(A | W, E) under observed treatment mechanism
#'  g0. If not specified, its form will be \code{Anodes ~ WEnodes}. If there are more than one expsosure, it fits a joint probability. See section 
#'  "Modeling P(A | W, E) for covariates (A, W, E)".
#' @param hform.gstar Character vector of regression formula for estimating the conditional density P(A | W, E) under user-supplied interventions  
#'  f_gstar1 or f_gstar2. If not specified, it use the same regression formula as used in hform.g0. 
#' @param lbound Value between (0,1) for truncation of predicted P(A | W, E). Default to 0.005
#' @param h.g0_GenericModel An object of \code{\link{GenericModel}} \pkg{R6} class containing the previously fitted models for P(A | W, E) under   
#'  observed treatment mechanism g0, one of the returns of \code{tmleCommunity} function. If known, predictions for P(A=a | W=w, E=e) under g0 are  
#'  based on the fitted models in \code{h.g0_GenericModel}.
#' @param h.gstar_GenericModel An object of \code{\link{GenericModel}} \pkg{R6} class containing the previously fitted models for P(A^* | W, E) under  
#'  intervention gstar, one of the returns of \code{tmleCommunity} function. If known, predictions for P(A=a | W=w, E=e) under gstar are based on 
#'  the fitted models in \code{h.gstar_GenericModel}.
#' @param TMLE.targetStep TMLE targeting step method, either "tmle.intercept" (Default) or "tmle.covariate". See "Details".
#' @param n_MCsims Number of simulations for Monte-Carlo analysis. Each simulation generates new exposures under f_gstar1 or f_gstar2 (if specified)  
#'  or f_g0 (if specified), with a sample size of nrow(data). Then these generated expsosures are used when fitting the conditional densities P(A | W, E).
#'  and estimating for \strong{IPTW} and \strong{GCOMP} under intervention f_gstar1 or f_gstar2. Note that deterministic intervention only needs one 
#'  simulation and stochastic intervention could use more simulation times such as 10 (Default to 1). 
#' @param CI_alpha Significance level (alpha) used in constructing a confidence interval. Default to 0.05.
#' @param rndseed Random seed for controlling sampling A under f_gstar1 or f_gstar2 (for reproducibility of Monte-Carlo simulations)
#' @param verbose Flag. If \code{TRUE}, print status messages. Default to \code{getOption("tmleCommunity.verbose")} (Default to \code{FALSE}). It 
#'  can be turned on by setting \code{options(tmleCommunity.verbose = TRUE)}.
#' 
#' @details
#' The estimates returned by \code{tmleCommunity} are of a treatment-specific mean, \eqn{E[Y_{g^*}]}, the expected community-level outcome if all 
#'  communities in the target population received the exposures generated under the user-supplied (deterministic or stochastic) intervention \eqn{g^*}. 
#'  
#' \code{data} must be a data frame (no matrix accepted) and doesn't support factor values (considering removing or recording such columns as strings), 
#'  \code{data} includes the following (optional) columns:
#'  \itemize{
#'   \item community-level and individual-level baseline covariate columns (\code{WEnodes}): can be any numeric data. Notice that W represents 
#'     individual-level covariates and E represent community-level covariates. 
#'   \item exposure columns (\code{Anodes}): can have more than one exposure and each exposure could be can be either binary, categorical or continuous. 
#'   \item outcome column (\code{Ynode}): can be any numeric data. If \code{Ynode} values are continuous, they may be automatically scaled. 
#'     See details for \code{Qbounds} below.  
#'   \item deterministic \code{Ynode} indicator column (\code{YnodeDet}): (optional) column that must be logical or binary. Only non-determinstic  
#'     \code{Ynode} values will be used in the final estimation step (e.g., \code{IPTW[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]}). 
#'   \item community identifier variable column (\code{communityID}): (optional) column that stores community identifier for hierarchical data 
#'     (or subject identifier if multiple observations for the same individual). Integer or character recommended (No factor allowed). 
#'     \code{NULL} means no hierarchical structure and all distinct individuals.
#'  }
#' 
#' \code{community.step} specifies how to read the structure of the input data (as hierarchical or non-hierarchical) and how to analyze the data. 
#'  community_level TMLE (\code{"community_level"}) is exactly analogous to the TMLE of the treatment specific individual-level mean outcome 
#'  (\code{"NoCommunity"}) with the trivial modification that the community rather than the individual serves as the unit of analysis.  
#'  \code{communityID} is needed when using \code{"community_level"}, \code{"individual_level"} and \code{"perCommunity"}. Lack of  
#'  \code{communityID} forces the algorithm to automatically pool data over all communities and treat it as non-hierarchical dataset (so force   
#'  \code{community.step} = \code{"NoCommunity"}). If \code{community.step} = \code{"individual_level"}, it incorporates with working models that 
#'  assume that each individual's outcome is known not to be affected by the covariates of other individuals in the same community (i.e., "no   
#'  covariate interference"). This strong assumption can be relaxed by integrating knowledge of the dependence structure among individuals  
#'  within communities (i.e., "weak covariate interference"). But currrently only the "no covariate interference" assumption is implemented. 
#'  If \code{community.step} = \code{"perCommunity"}), run a single TMLE on each community and calculate a (weighted) mean outcome for the J  
#'  communities.
#' 
#' \code{pooled.Q} is in regard to incorporate the working model of "no covariate inference" in community-level \code{TMLE} (and the 
#'  corresponding \code{IPTW} and \code{GCOMP}) although the working model is not assumed to hold. In other words, when \code{community.step} 
#'  = \code{"community_level"}, if \code{pooled.Q} = \code{TRUE}, add pooled individual-level regressions as candidates in the Super Learner  
#'  library for initial estimation of the outcome mechanism. If \code{pooled.Q} = \code{FALSE}, both outcome and treatment mechanisms are
#'  estimated on the community-level (no use of individual-level information). 
#'  
#' \code{Qform} should be \code{NULL}, in which cases all parent nodes of Y node will be used as regressors, or a character vector that can  
#'  be coerced to class \code{"formula"}. If \code{Qestimator} (an argument in \code{tmleCom_Options}) is \code{"speedglm__glm"} (or  
#'  \code{"speedglm__glm"}), then \code{speedglm} (or \code{glm}) will be called using the components of \code{Qform}. If \code{Qestimator} 
#'  is \code{"SuperLearner"}, then \code{SuperLearner} will be called after a data frame is created using \code{Qform}, based on the specified
#'  algorithms in \code{SL.library} (an argument in \code{tmleCom_Options}); If \code{Qestimator} is \code{"h2o__ensemble"}, then \code{h2o} and 
#'  \code{h2oEnsemble} will be called after a H2OFrame dataset is creating using \code{Qform}, based on specified algorithms in \code{h2olearner}
#'  and \code{h2ometalearner}. See "Arguments" in \code{tmleCom_Options}.
#' 
#' \code{hform.g0} and \code{hform.gstar} should also be \code{NULL}, in which cases all parent nodes of A node(s) will be used as regressors, or 
#'  a character vector that can be coerced to class \code{"formula"}. It follows the same rules applied to \code{Qform} except it's now based on
#'  \code{gestimator} (an argument in \code{tmleCom_Options}). See "Arguments" in \code{tmleCom_Options}.
#' 
#' \code{Qbounds} can be used to bound continuous outcomes Y. If \code{Qbounds} not specified (\code{NULL}), it will be set to the range of Y, 
#'  widened by 10\% of the minimum and maximum. That is, [0.9*\code{min(Y)}, 1.1*\code{max(Y)}]. If specified, then Y will be truncated at the min 
#'  max values of \code{Qbounds}, and then scaled to be in [0, 1] by (Y - \code{min(Qbound)})/(\code{diff(Qbound)}). Statistical inferences  
#'  and for the transformed outcome can be directly translated back into statistical inference for the unscaled outcome. Once \code{Qbounds} finish
#'  bounding the observed outcomes, it will be set to (1 - \code{alpha}, \code{alpha}) and used to bound the predicted values for the initial outcome 
#'  mechanism. Thus, \code{alpha} needs to be between (0, 1), otherwise reset to 0.995. Besides, \code{lbound} can be used to truncate the weights 
#'  \code{h_gstar/h_gN}, that is, [0, 1/\code{lbound}].
#' 
#' \code{TMLE.targetStep} specifies how to use weights \code{h_gstar/h_gN} in the \strong{TMLE} targeting step. If \code{tmle.intercept} (default),   
#'  it performs the weighted intercept-based TMLE that runs a intercept-only weighted logistic regression using offsets \code{logit(Qstar)} and  
#'  weights \code{h_gstar/h_gN} and so no covariate. If \code{tmle.covariate}, it performs the unweighted covariate-based TMLE that run a 
#'  unweighted logistic regression using offsets \code{logit(Qstar)} and a clever covariate \code{h_gstar/h_gN}. 
#' 
#' @section IPTW estimator:
#' The IPTW estimator is based on the TMLE weights \eqn{h_{g^*}(A^*,W,E)/h_{g}(A,W,E)}, which is equivalent to \eqn{P_{g^*}(A^*|W,E)/P_{g}(A|W,E)} 
#'  and is defined as the the weighted average of the observed outcomes Y. The following algorithm shows a general template of the community-level 
#'  \code{IPTW}: 
#'  \itemize{
#'   \item As described in the following section ("Modeling \code{P(A | W, E)} for covariates \code{(A, W, E)}"), the first step is to construct an 
#'     estimator \eqn{P_{\hat{g}^{c}}(A | W, E)} of the density for the common (in j) conditional distribution of \eqn{A} given \eqn{W, E}, that is 
#'     \eqn{P_{g_0^{c}}(A | W, E)} for common (in j) community-level covariates (A, W, E). 
#'   \item Implementing the same modeling & fitting algorithm to construct an estimator \eqn{P_{\hat{g}^{c*}}(A^* | W, E)} of the density for  
#'     the common (in j) conditional distribution of \eqn{A^*} given \eqn{(W, E)}, that is \eqn{P_{g_0^{c*}}(A^* | W, E)} for common (in j)
#'     community-level covariates (\eqn{A^*, W, E}) where \eqn{A^*} is determined by the user-supplied stochastic intervention \code{f_gstar1} or 
#'     \code{f_gstar2}, given the observed baseline covariates \eqn{(W, E)}. 
#'   \item Given observed J independent communities \eqn{\textbf{O}_j = (E_j, \textbf{W}_j, A_j, Y^c_j: j = 1, ..., J)}, the IPTW estimator is given by:
#'     \deqn{\psi^{I}_{IPTW, n}=\frac{1}{J}\sum_{j=1}^{J}Y^c_j\frac{P_{\hat{g}^{c*}}(A^*_j|\textbf{W}_j,E_j)}{P_{\hat{g}^{c*}}(A_j|\textbf{W}_j,E_j)}}
#'  }
#'
#' For individual-level IPTW, it reads the input data as \eqn{O_{i,j} = (E_j, W_{i,j}, A_j, Y_{i,j}: j = 1, ..., J; i = 1, ..., n_j)} and incorporates 
#'  working model that assumes no covariate interference, weighing each individual within one community by \eqn{\alpha_{i,j}}, where the IPTW estimator
#'  is given by: \deqn{\psi^{II}_{IPTW, n}=\frac{1}{J}\sum_{j=1}^{J}\sum_{i=1}^{n_j}\alpha_{i,j}Y_{i,j}\frac{P_{g^{*}}(A^*_j|W_{i,j}, E_j)}
#'   {P_{g^{*}}(A_j | W_{i,j}, E_j)}}
#' 
#' @section TMLE estimator:
#' The TMLE estimator is based on the updated model prediction \eqn{\bar{Q}^*(A, W, E)} and is defined by the G-formula. The following algorithm  
#'  shows a general template of the community-level \code{TMLE}:
#'  \itemize{ 
#'   \item The first step is exactly the same as \code{IPTW}: construct two density estimators and use the ratio of them as the weights 
#'     \eqn{P_{g^*}(A^*|W,E)/P_{g}(A|W,E)} in the targeting step. 
#'   \item Construct an initial estimator \eqn{\hat{\bar{Q}}^c(A | W, E)} of the common (in j) conditional distribution of \eqn{Y^c} given 
#'     \eqn{(A, W, E)} and update \eqn{\hat{\bar{Q}}^{c*}(A | W, E)} for \eqn{\hat{\bar{Q}}^c(A | W, E)} by weights calculated in the first step.
#'   \item The TMLE estimator is defined as the following substitution estimator:
#'     \deqn{\psi^{I}_{TMLE, n}=\frac{1}{J}\sum_{j=1}^{J}\int_{a}\hat{\bar{Q}}^{c*}(a, \textbf{W}_j, E_j)g^{c*}(a|\textbf{W}_j, E_j)d\mu(a)}
#'  }
#'
#' For individual-level TMLE, its estimator is obtained as:
#'  \deqn{\psi^{II}_{TMLE,n}=\frac{1}{J}\sum_{j=1}^{J}\sum_{i=1}^{n_j}\alpha_{i,j}\int_{a}\hat{\bar{Q}}^*(a, W_{i,j}, E_j)g^*(a|W_{i,j}, E_j)d\mu(a)}
#' 
#' @section GCOMP estimator:
#' The GCOMP estimator is similar to the the TMLE estimator except it uses the untargeted (initial) model \eqn{\hat{\bar{Q}}^c(A|W,E)} instead of 
#'  its targeted version \eqn{\hat{\bar{Q}}^{c*}(A | W, E)}, for example the community-level GCOMP estimator is given by: 
#'  \deqn{\psi^{I}_{GCOMP,n}=\frac{1}{J}\sum_{j=1}^{J}\int_{a}\hat{\bar{Q}}^{c}(a, \textbf{W}_j, E_j)g^{c*}(a|\textbf{W}_j, E_j)d\mu(a)}
#  Similary the individual estimator can be written as:
#  \deqn{\psi^{II}_{GCOMP,n}=\frac{1}{J}\sum_{j=1}^{J}\sum_{i=1}^{n_j}\alpha_{i,j}\int_{a}\hat{\bar{Q}}(a, w_{i,j}, E_j)g^*(a|w_{i,j}, E_j)d\mu(a)}
#'
#' @section Modeling \eqn{P(A | W, E)} for covariates \eqn{(A, W, E)}:
#' For simplicity (and without loss of generality), we now suppose that there is no hierarchical structure in data and are interested in finding 
#'  an non-parametric estimator of the common (in i) \strong{individual-level} exposure mechanism \eqn{g_0(A|W)}, or the commom multivariate 
#'  joint conditional probability model \eqn{P_{g_0}(A|W)}, where the exposures and baseline covariates (\eqn{(A,W)=(A_i, W_i: i=1,...,n)}) denote 
#'  the random variables drawn jointly from distribution \eqn{H_0(A, W)} with denisty \eqn{h_0(a, w) \equiv g_0(a|w)q_{W,0}(w)} and \eqn{q_{W,0}(W)}
#'  denotes the marginal density of the baseline covariates \eqn{W} specified by the regression formula \code{hform.g0}
#'  (Notice that an non-parametric estimator of the model \eqn{P_{g^*_0}(A|W)}) is similar, except that now the exposures and baseline covariates 
#'  (\eqn{(A^*,W)=(A^*_i, W_i: i=1,...,n)}) are randomly drawn from \eqn{H^*_0(A, W)} with density \eqn{h^*_0(a, w) \equiv g^*_0(a|w)q_{W,0}(w)}, 
#'  where \eqn{A^*} is determined by the user-supplied (stochastic) intervention \code{f_gstar1} or \code{f_gstar2} and \eqn{q_{W,0}(W)} denotes
#'  the marginal density of the baseline covariates \eqn{W} specified by the regression formula \code{hform.gstar}. Thus, the fitting algorithm
#'  for \eqn{P_{g^*_0}(A|W)}) is equivalent for \eqn{P_{g_0}(A|W)}).
#' 
#' Note that \eqn{A} can be multivariate (i.e., \eqn{(A(1), ..., A(M))}) and each of its components \eqn{A(m)} can be either binary, categorical  
#'  or continuous. The joint probability model for \eqn{P(A|W)=P(A(1),...,A(M)|W)} can be factorized as a sequence \eqn{P(A(1)|W) \times 
#'  P(A(2)|W,A(1)) \times ... \times P(A(M)|W, A(1),...,A(M-1))}, where each of these conditional probability models \eqn{P(A(m)|W,A(1),...,A(m-1))}
#'  is fit separately, depending on the type of the m-specific outcome variable \eqn{A(m)}. For binary \eqn{A(m)}, the conditional probability 
#'  \eqn{P(A(m)|W,A(1),...,A(m-1))} will be esimtated by a user-specific library of candidate algorithms, including parametric estimators such as 
#'  logistic model with only main terms, and data-adaptive estimator such as super-learner algorithms. For continuous (or categorical) \eqn{A(m)},
#'  consider a sequence of values \eqn{\delta_1, \delta_2,...,\delta_{K+1}} that span the range of \eqn{A} and define \eqn{K} bins and the 
#'  corresponding \eqn{K} bin indicators (\eqn{B_1,...,B_K}), in which case each bin indicator \eqn{B_k \equiv [\delta_k, \delta_{k+1})} is used 
#'  as an binary outcome in a seperate user-specific library of candidate algorithms, with predictors given by \eqn{(W, A(1),..., A(m-1))}. That is  
#'  how the joint probability \eqn{P(A|W)} is factorized into such an entire tree of binary regression models.
#'
#' For simplicity (and without loss of generality), we now suppose \eqn{A} is univariate (i.e., M=1) and continuous and a general template of an 
#'  fitting algorithm for \eqn{P_{g_0}(A|W)} is summrized below:
#'
#' \enumerate{
#' \item Consider the usual setting in which we observe n independently and identically distributed copies \eqn{o_i=(w_i, a_i, y_i: i=1,..,n)} of
#'   the random variable \eqn{O=(W, A, Y)}, where the observed \eqn{(a_i: i=1...,n)} are continuous. 
#'
#' \item As described above, consider a sequence of \eqn{K+1} values that span the support of \eqn{A} values into \eqn{K} bin intervals \eqn{\Delta} 
#'   = (\eqn{\delta_1, \delta_2,...,\delta_{K+1}}) so that any observed data point \eqn{a_i} belongs to one interval within R, in other words, 
#'   for each possible value \eqn{a \in A} (even if it's not in the observed \eqn{(a_i:i)}), there always exists a \eqn{k \in {1, ...,K}} such  
#'   that \eqn{\delta_{k}\leq a<\delta_{k+1}}, and the length (bandwidth) of the interval can be defined as \eqn{bw_{k}=[\delta_{k+1}-\delta_{k})}. 
#'   Then let the mapping \eqn{S(a)\in \{1,2,..,K\}} denote a unique index of the indicator in \eqn{\Lambda} that \eqn{a} falls in, where \eqn{S(a)=k}
#'   if \eqn{a\in [\delta_{k},\delta_{k+1})}, namely, \eqn{\delta_{S(a)} \leq a < \delta_{S(a)+1}}. Moreover, we use \eqn{b_k} to denote a binary 
#'   indicator of whether the observed \eqn{a} belongs to bin \eqn{k} (i.e., \eqn{b_k\equiv I(S(a)=k)} for all \eqn{k\leq S(a)}; \eqn{b_k\equiv}
#'   \code{NA} for all \eqn{k>S(a))}. This is similar to methods for censored longitudinal data, which code observations as \code{NA} (censored or 
#'   missing) once the indicator \eqn{b_k} jumps from 0 to 1. Since \eqn{a} is a realization of the random variable \eqn{A} for one individual, 
#'   the corresponding random binary indicators of whether \eqn{A} belongs to bin \eqn{k} can be denoted by \eqn{B_k:k=1,..,=K} where 
#'   \eqn{B_k \equiv I(S(A)=k)} for all \eqn{k\leq S(A)}; \eqn{B_k\equiv}\code{NA} for all \eqn{k>S(A)}.
#'
#' \item Then for each k = 1,...,K, a binary nonparametric regression is used to estimate the conditional probability \eqn{P(B_k=1|B_{k-1}=0,W)}, 
#'   which corresponds to the probability of \eqn{B_k} jumping from 0 to 1, given \eqn{B_{k-1}=0} and tbe baseline covariates \eqn{W}. Note tha 
#'   for each k, the corresponding nonparametric regression model is fit only among observations that are uncensored (i.e., still at risk of  
#'   getting \eqn{B_{k}=1} with \eqn{B_{k-1}=0}). Note the above conditional probability \eqn{P(B_k=1|B_{k-1}=0,W)} is equivalent to 
#'   \eqn{P(A\in [\delta_{k}, \delta_{k+1}) | A\geq \delta_{k+1}, W)}, which is the probability of \eqn{A} belongs to the interval 
#'   \eqn{[\delta_{k},\delta_{k+1})}, conditional on \eqn{A} doesn't belong to any intervals before \eqn{[\delta_{k}, \delta_{k+1})} and \eqn{W}.
#'   Then the discrete conditional hazard function for each k is defined as a normalization of the conditional probability using the 
#'   corresponding interval bandwidth \eqn{bw_{k}}: 
#'   \eqn{\lambda_k(A,W)=\frac{P(B_k=1|B_{k-1}=0,W)}{bw_k}=\frac{P(A\in[\delta_{k},\delta_{k+1})|A\geq \delta_{k+1},W)}{bw_k}}
#'
#' \item Finally, for any given observation \code{(a,w)}, we first find out the interval index \eqn{k} to which \eqn{a} belongs (i.e., 
#'   \eqn{k=S(a)\in{1,...,K}}). Then the discretized conditional density of \eqn{P(A=a|W=w)} can be evaluated by 
#'   \eqn{\lambda_k(A, W){\times}[\prod_{j=1}^{k-1}(1-\lambda_j(A, W))]}, which corresponds to the conditional probability of \eqn{a} 
#'   belongs to the interval \eqn{[\delta_{k},\delta_{k+1})} and doesn't belong to any intervals before, given \eqn{W}.
#' }
#'
#' @section Three methods for choosing bin (interval) locations for a univariate and continuous variable \eqn{A}:
#' Note that the choice of the values \eqn{\delta_k (k=1,..,K)} implies defining the number and positions of the bins. First, a cross- 
#'  validation selector can be applied to data-adaptively select the candidate number of bins, which minimizes variance and maximizes precision 
#'  (Don't recommend too many bins due to easily violating the positivity assumption). Then, we need to choose the most convenient locations 
#'  (cuttoffs) for the bins (for fixed K). There are 3 alternative methods that use the histogram as a graphical descriptive tool to define 
#'  the bin cutoffs \eqn{\Delta=(\delta_1,...,\delta_K,\delta_{K+1})} for a continuous variable \code{A}. In \pkg{tmleCommunity} package, the 
#'  choice of methods \code{bin.method} together with the other discretization arguments in function \code{tmleCom_Options()} such as \code{nbins}  
#'  (total number of bins) and \code{maxNperBin} (the maximum number of observations in each bin), can be used to define the values of bin  
#'  cutoffs. See \code{help(tmleCom_Options)} for more argument details. 
#' 
#' \itemize{
#'  \item \code{equal.mass}: The default discretization method, equal mass (aka equal area) interval method, set by passing an argument
#'    \code{bin.method="equal.mass"} to \code{tmleCom_Options()} prior to calling the main function \code{tmleCommunity()}. The interval are 
#'    defined by spanning the support of \eqn{A} into non-equal length of bins, each containing (approximately) the same number of observations. 
#'    It's data-adaptive since it tends to be wide where the population density is small, and narrow where the density is large. If \code{nbins} 
#'    is \code{NA} (or is smaller than \code{n/maxNperBin}), \code{nbins} will be (re)set to the interger value of \code{n/maxNperBin} where 
#'    \code{n} is the total number of observations in \eqn{A}, and the default setting of \code{maxNperBin} is 500 observations per interval.
#'    This method could identify spikes in the density, but oversmooths in the tails and so couldn't discover outliers. 
#'    
#'  \item \code{equal.len}: equal length interval method, set by passing an argument \code{bin.method="equal.len"} to 
#'    \code{tmleCom_Options()} prior to calling \code{tmleCommunity()}. The intervals are defined by spanning the support of \eqn{A} into 
#'    \code{nbins} number of equal length of bins. This method describes the tails of the density and identifies outliers well, but  
#'    oversmooths in regions of high density and so is poor at identifying sharp peaks.
#'
#'  \item \code{dhist}: (named for diagonally cut histogram) Combination of equal-area equal length and equal mass method, set by passing an 
#'    argument \code{bin.method="dhist"} to \code{tmleCom_Options()} prior to calling \code{tmleCommunity()}. For consistency, We choose the 
#'    slope a = \eqn{5\times \textrm{IQR}(A)} as suggested by Denby and Mallows ("Variations on the Histogram" (2009)). For more details,   
#'    please also see this paper.
#' }
#'
#' @return \code{tmleCommunity} returns an object of class "\code{tmleCommunity}", which is a named list containing the estimation results
#'  stored by the following 3 elements:
#'  \itemize{
#'  \item \code{EY_gstar1} - a list with estimates of the mean counterfactual outcome under (deterministic or stochastic) intervention  
#'    function \code{f_gstar1} (\eqn{E[Y_{g^*_1}]}).
#'  \item \code{EY_gstar2} - a list with estimates of the mean counterfactual outcome under (deterministic or stochastic) intervention   
#'    function \code{f_gstar2} \eqn{(E[Y_{g^*_2}])}; or \code{NULL} if \code{f_gstar2} not specified.
#'  \item \code{ATE} - a list with estimates of additive treatment effect (\eqn{E[Y_{g^*_1}]} - \eqn{E[Y_{g^*_2}]}) under two interventions
#'     functions \code{f_gstar1} VS \code{f_gstar2}; or \code{NULL} if \code{f_gstar2} not specified.
#' }
#' Each element in the returned \code{tmleCommunity} object is itself a list containing the following 8 items:
#'  \itemize{
#'  \item \code{estimates} - matrix, 3\eqn{\times}1, storing 3 algorithm estimates of the target parameter (population community-level counterfactual 
#'    mean under (deterministic or stochastic) intervention), including \code{TMLE}, \code{IPTW} and \code{GCOMP}.
#'  \item \code{vars} - matrix, 3\eqn{\times}1, storing 3 influence-curve based asymptotic variance estimates for \code{TMLE}, \code{IPTW} and 
#'    \code{GCOMP}. Notice, all IC-based statistical inference for \strong{GCOMP} is not accurate (Just for reference). See explanation in \code{IC}.
#'  \item \code{CIs} - matrix, 3\eqn{\times}2, storing 3 confidence interval estimates at \code{CI_alpha} level, for \code{TMLE}, \code{IPTW} and 
#'    \code{GCOMP}. The first column contains the lower bounds and the second column contains the upper bounds.
#'  \item \code{tstat} - matrix, 3\eqn{\times}1, storing 3 test statistics.
#'  \item \code{pval} - matrix, 3\eqn{\times}1, storing 3 p-values. It's based on the Student's T distribution if the number of communities  
#'    (or the number of individuals if no hierarchical structure) is less than 41, otherwise based on the Z normal distribution.
#'  \item \code{IC} - data frame, nobs\eqn{\times}3, the first column contains the influence curves (ICs) for \code{TMLE} estimate, the second column
#'    contains the ICs for \code{IPTW} estimate, and the third column contains the ICs for \code{GCOMP} estimate (not accurate since it's based on ICs
#'    for \code{TMLE} estimate without updating step).  
#'  \item \code{h.g0_GenericModel} - An object of \code{\link{GenericModel}} \pkg{R6} class, storing the model fits for P(A | W, E) under  
#'    observed exposure mechanism \code{g0}. This can be used in \code{\link{tmleCommunity}} (See \strong{Arguments}).
#'  \item \code{h.gstar_GenericModel} - An object of \code{\link{GenericModel}} \pkg{R6} class, storing the model fits for P(A | W, E)  
#'    under intervention \code{f_gstar1} or \code{f_gstar2}. This can be used in \code{\link{tmleCommunity}} (See \strong{Arguments}).
#' }
#' 
#' Estimations are based on either community-level or individual-level analysis. Each analysis currently implements 3 estimators:
#'  \itemize{
#'  \item \code{tmle} - Either weighted intercept-based TMLE based on weights \code{h_gstar/h_gN} (default choice)
#'    or unweighted covariate-based TMLE based on a covariate \code{h_gstar/h_gN}.  
#'  \item \code{iptw} - IPTW (Horvitz-Thompson) estimator based on the \code{TMLE} weights h_gstar/h_gN.
#'  \item \code{gcomp} - Maximum likelihood based G-computation substitution estimator.
#' }
#' @seealso \code{\link{tmleCommunity-package}} for the general description of the package,
#'  
#'  \code{\link{tmleCom_Options}} for additional parameters that control the estimation in \code{tmleCommunity},
#' 
#'  \code{\link{DatKeepClass}} for details about storing, managing, subsetting and manipulating the input data,
#' 
#'  \code{\link{indSample.iid.cA.cY_list}} for an example of a continuous exposure and its estimation,
#' 
#'  \code{\link{BinaryOutModel}}, \code{\link{RegressionClass}}, \code{\link{GenericModel}}, \code{\link{MonteCarloSimClass}}
#' 
#' @example tests/examples/6_tmleCommunity_examples.R
#' @export
tmleCommunity <- function(data, Ynode, Anodes, WEnodes, YnodeDet = NULL, obs.wts = c("equal.within.pop", "equal.within.community"), 
                          community.step = c("NoCommunity", "community_level", "individual_level", "perCommunity"),
                          communityID = NULL, community.wts = c("size.community", "equal.community"), pooled.Q = FALSE, 
                          f_g0 = NULL, f_gstar1, f_gstar2 = NULL, Qform = NULL, Qbounds = NULL, alpha = 0.995,                                                      
                          fluctuation = "logistic", hform.g0 = NULL, hform.gstar = NULL, lbound = 0.005, 
                          h.g0_GenericModel = NULL, h.gstar_GenericModel = NULL, 
                          TMLE.targetStep = c("tmle.intercept", "tmle.covariate"),
                          n_MCsims = 1, CI_alpha = 0.05, rndseed = NULL, verbose = getOption("tmleCommunity.verbose")) {
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  oldverboseopt <- getOption("tmleCommunity.verbose")
  options(tmlenet.verbose = verbose)
  gvars$verbose <- verbose
  if (verbose) { message("Running tmleCommunity with the following settings from tmleCom_Options(): "); str(gvars$opts) }
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (missing(Ynode)) Ynode <- NULL
  savetime.fit.hbars <- getopt("savetime.fit.hbars")  # One of conditions for skip g0 & gstar fitting procedure and directly set h_gstar_h_gN = 1 
  if (is.character(obs.wts)) {
    obs.wts <- obs.wts[1]
    if (obs.wts == "equal.within.pop") { # weigh individuals in the entire dataset equally so big community gets bigger total weight
      obs.wts <- rep(1, NROW(data))
    } else if (obs.wts == "equal.within.community") { # weigh individuals in each community equally and weigh communities equally
      obs.wts <- rep(as.vector(1/table(data[, communityID])), as.vector(table(data[, communityID])))
    }
  }
  if (is.character(community.step)) {
    community.step <- community.step[1]
    if ((community.step != "NoCommunity") && !is.data.frame(community.wts) && !is.null(communityID)) {
      community.wts <- community.wts[1]
      community.wts.df <- as.data.frame(matrix(0L, nrow = length(unique(data[, communityID])), ncol = 2))
      colnames(community.wts.df) <- c("id", "weights")
      community.wts.df[, 1] <- names(table(data[, communityID]))
      if (community.wts == "equal.community") { # weigh each community equally
        community.wts.df[, 2]  <- rep(1, length(unique(data[, communityID]))) 
      } else if (community.wts == "size.community") { # weigh each community by its number of obs - larger community has larger weight
        community.wts.df[, 2] <- as.vector(table(data[, communityID]))
      }
      community.wts <- community.wts.df; community.wts.df <- NULL
    } else if (community.step == "NoCommunity" || is.null(communityID)) { # Don't need community.wts anymore, just need to pass CheckInputs()
       community.wts <- NULL
    }
  } 
  if (!is.null(Qform) && !is.null(Ynode)) {
    Qform <- paste(Ynode, substring(Qform, first = as.numeric(gregexpr("~", Qform))))
    if (verbose) {
      message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+% Ynode)
      message("Thus the Qform becomes " %+% Qform)
    }
  }
  if (!is.null(Qform) && is.null(Ynode)) {
    Ynode <- LhsVars(Qform)[1]
    message("Setting the Ynode to: " %+% Ynode)
  }
  TMLE.targetStep <- TMLE.targetStep[1]
  
  ## Check if any unexpected inputs
  if (!(TMLE.targetStep %in% c("tmle.intercept", "tmle.covariate"))) 
    stop("TMLE.targetStep argument must be either 'tmle.intercept' or 'tmle.covariate'")
  if (!(community.step %in% c("NoCommunity", "community_level", "individual_level", "perCommunity"))) 
    stop("community.step argument must be one of 'NoCommunity', 'community_level', 'individual_level' and 'perCommunity'")
  if ((community.step %in% c("community_level", "individual_level", "perCommunity")) & is.null(communityID)) {
    messageMSg <- c("'communityID' is needed when using 'community_level', 'individual_level' and 'perCommunity'. Lack of 'communityID'",
                    "forces the algorithm to automatically pool data over all communities and treat it as non-hierarchical dataset. \n",
                    "In other words, we simply treat community.step = 'NoCommunity' and pooled.Q = FALSE")
    message(messageMSg[1] %+% messageMSg[2] %+% messageMSg[3])
    community.step <- "NoCommunity"; pooled.Q <- FALSE
  } 
  if (!(community.wts %in% c("equal.community", "size.community")) && !is.data.frame(community.wts) && (community.step != "NoCommunity")) {
    stop("Currently only numeric values, 'equal.community' and 'size.community' are supported for community.wts")
  }
  nodes <- list(Ynode = Ynode, Anodes = Anodes, WEnodes = WEnodes, communityID = communityID)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, Qform, hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, community.wts, f_gstar1, f_gstar2)) stop()
  # ensure that the column names are defined as "id" and "weights"
  if (community.step != "NoCommunity") colnames(community.wts) <- c("id", "weights")  
  maptoYstar <- fluctuation=="logistic"  # if TRUE, cont Y values shifted & scaled to fall b/t (0,1)
  
  #----------------------------------------------------------------------------------
  # DEFINING (OPTIONAL) REGRESSION FORMS 
  #----------------------------------------------------------------------------------
  Q.sVars <- define_regform(as.formula(Qform), Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "WEnodes")])
  h.g0.sVars <- define_regform(as.formula(hform.g0), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes["WEnodes"])
  if (!is.null(hform.gstar)) {
    h.gstar.sVars <- define_regform(as.formula(hform.gstar), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes["WEnodes"])
  } else {
    h.gstar.sVars <- h.g0.sVars
  }
  
  if (verbose) {
    print("Input regression Qform (E(Y|A,W,E)): " %+% Qform)
    print("Derived regression Qform (E(Y|A,W,E)):"); str(Q.sVars)
    print("Input regression hform.g0 (P(A|W,E) under g0): " %+% hform.g0)
    print("Derived regression hform.g0 (P(A|W,E) under g0): "); str(h.g0.sVars)
    print("Input regression hform.gstar (P(A|W,E) under g.star): " %+% hform.gstar)
    print("Derived regression hform.gstar (P(A|W,E) under g.star): "); str(h.gstar.sVars)
  }
  
  ## Create data based on community.step, then based on Qform, hform.g0 and hform.gstar, in case of interaction or higher-order term.
  if (community.step == "community_level" && !pooled.Q) { 
    f_gstar1.agg <- f_gstar2.agg <- FALSE
    # if running entire TMLE algorithm at cluster-level without a pooled individual-level regression on outcome, aggregate data now 
    # Also need to aggregate f_gstar1 and f_gstar2 if they are vectors of length NROW(data) 
    if (is.VecMatDf(f_gstar1) && NROW(as.data.frame(f_gstar1)) == NROW(data)) { f_gstar1.agg <- TRUE; data <- cbind(data, f_gstar1) }
    if (is.VecMatDf(f_gstar2) && NROW(as.data.frame(f_gstar2)) == NROW(data)) { f_gstar2.agg <- TRUE; data <- cbind(data, f_gstar2) }
    data <- aggregate(x = data, by=list(newid = data[, communityID]), mean) # [, 2 : (ncol(data)+1)] # Don't keep the extra ID column
    colname.allNA <- colnames(data)[colSums(is.na(data)) == NROW(data)]  # columns with all NAs after aggregation, due to non-numeric values
    if (length(colname.allNA) != 0) {
      data <- data[, -which(colnames(data) %in% colname.allNA)] # Remove columns that contain only NAs inside
      names(data)[which(colnames(data) == "newid")] <- communityID  # change 'newid' back to the communityID name
      warning(paste(colname.allNA, collapse = ', ') %+% " is(are) removed from the aggregated data due to all NAs in the column(s).")
      warning("Suggestion: convert the non-numeric values to numeric, e.g., create dummy variables for each category/ string/ factor.")
    }
    # ensure that weights match with their corresponding communities
    obs.wts <- community.wts[match(data[, communityID], community.wts[, "id"]), "weights"]  
    if (f_gstar1.agg) { if (is.vector(f_gstar1)) { f_gstar1 <- data$f_gstar1 } else {f_gstar1 <- data[, names(f_gstar1)]} } 
    if (f_gstar2.agg) { if (is.vector(f_gstar2)) { f_gstar2 <- data$f_gstar2 } else {f_gstar2 <- data[, names(f_gstar2)]} } 
  }
  
  if (!is.null(c(Qform, hform.g0, hform.gstar))) {
    allcovRHS <- unique(unlist(lapply(c(Qform, hform.g0, hform.gstar), FUN = function(x) { strsplit(deparse(as.formula(x)[[3]]), " \\+ ")[[1]] })))
    merged.form <- reformulate(allcovRHS, response = NULL)  # Reformulate a formula including all legitimate character of the RHS in 3 formulae
    if (any(!allcovRHS %in% unique(c(unlist(nodes), names(data))))) {
      ExtraDat <- as.data.frame(model.matrix(merged.form, data = data))
      data <- cbind(data, ExtraDat[, setdiff(names(ExtraDat), c(names(data), "(Intercept)")), drop = FALSE])
      ExtraDat <- NULL
    }
  } else {
    merged.form <- NULL
  }
  nodes <- append(nodes, list(Crossnodes = setdiff(names(data), Reduce(c, nodes))))
  # Why to keep variables that are not indicated in the node list (i.e. Crossnodes): when creating A^* under g.star (delta function)  
  # it's possible to use variablas that are not used in Qform and gform. 
  
  if (community.step %in% c("NoCommunity", "community_level", "individual_level")) {
    ## Create an R6 object that stores and manages the input data, later passed on to estimation algorithm(s)
    inputYs <- CreateInputs(data[, Ynode], Qbounds, alpha, maptoYstar)
    data[, Ynode] <- inputYs$Ystar
    OData.ObsP0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
    OData.ObsP0$addYnode(YnodeVals = inputYs$Ystar)
    OData.ObsP0$addObsWeights(obs.wts = obs.wts)
    nobs <- OData.ObsP0$nobs
    if (is.null(YnodeDet)) {
      determ.Q <- rep_len(FALSE, nobs)
    } else {
      determ.Q <- (data[, YnodeDet] == 1)
    }
    if (length(unique(obs.wts)) > 1 && any(unlist(OData.ObsP0$type.sVar[Anodes]) != "binary")) {
      if (verbose) warning("obs.wts are currently implemented on binary A. The results for non-binary A with weights may be unrealiable.")
    }
    
    #----------------------------------------------------------------------------------
    # Defining and estimating outcome mechanism E(Y|A, E, W)
    #----------------------------------------------------------------------------------
    if (verbose) {
      message("================================================================")
      message("fitting E(Y|A,W,E):= ", "P(" %+% nodes$Ynode %+% "=1 | " %+% paste(Q.sVars$predvars, collapse = ",") %+% ")")
      message("================================================================")
    }
    Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, 
                                predvars = Q.sVars$predvars, 
                                subset_vars = !determ.Q, 
                                estimator = getopt("Qestimator"))
    model.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData.ObsP0, savespace = TRUE)
    
    #----------------------------------------------------------------------------------
    # Create an list with model estimates, data & other information that is passed on to treatment estimation procedure
    #----------------------------------------------------------------------------------
    estinfo_list <- list(
      data = data, 
      nodes = nodes,
      communityID = communityID,
      community.step = community.step,
      pooled.Q = pooled.Q,
      TMLE.targetStep = TMLE.targetStep,
      obs.wts = obs.wts, 
      community.wts = community.wts,
      Qbounds = inputYs$Qbounds,
      lbound = lbound,
      merged.form = merged.form, 
      model.Q.init = model.Q.init,
      Q.sVars = Q.sVars,
      f.g0 = f_g0,
      h.g0.sVars = h.g0.sVars,
      h.gstar.sVars = h.gstar.sVars,
      h.g0_GenericModel = h.g0_GenericModel,
      h.gstar_GenericModel = h.gstar_GenericModel,
      savetime.fit.hbars = savetime.fit.hbars,
      n_MCsims = n_MCsims
    ) 
    estinfo_list_g1 <- append(estinfo_list, list(f.gstar = f_gstar1))
    if (!is.null(f_gstar2)) { estinfo_list_g2 <- append(estinfo_list, list(f.gstar = f_gstar2)) }
    
    #----------------------------------------------------------------------------------
    # Running MC evaluation for substitution TMLE estsimators
    #----------------------------------------------------------------------------------
    # Incl. estimate treatment mechanism f(a|E, W)) and clever covariates & targeting step
    tmle_gstar1_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g1)
    if (!is.null(f_gstar2)) {
      tmle_gstar2_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g2)
    } else {
      tmle_gstar2_out <- NULL
    }
    
    #----------------------------------------------------------------------------------
    # Create output list (estimates, as. variances, CIs)
    #----------------------------------------------------------------------------------
    EY_gstar1 <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, tmle_g_out = tmle_gstar1_out)
    EY_gstar2 <- NULL
    ATE <- NULL	
    if (!is.null(f_gstar2)) {
      EY_gstar2 <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, tmle_g_out = tmle_gstar2_out)
      ATE <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, tmle_g_out = tmle_gstar1_out, tmle_g2_out = tmle_gstar2_out)
    }
    if (verbose) {
      message("######################################################################################")
      message("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.")
      message("######################################################################################")
    }
  } else if (community.step == "perCommunity") {
    #----------------------------------------------------------------------------------
    # Create matrices to store substitution estsimators of all communities
    #----------------------------------------------------------------------------------
    communityList <- unique(data[, communityID])
    est.communities_gstar1 <- matrix(0L, nrow = length(communityList), ncol = 3)
    colnames(est.communities_gstar1) <- c("TMLE", "IPTW", "MLE")
    wts.communities_gstar1 <- fWi.communities_gstar1 <- matrix(0L, nrow = 0, ncol = 1)
    colnames(wts.communities_gstar1) <- c("h_wts")
    colnames(fWi.communities_gstar1) <- c("fWi_Qinit")
    QY.communities_gstar1  <- matrix(0L, nrow = 0, ncol = 2)
    colnames(QY.communities_gstar1) <- c("QY.init", "QY.star")
    obs.wts.communities <- c()
    f_gstar1.splt <- is.VecMatDf(f_gstar1) && NROW(as.data.frame(f_gstar1)) == NROW(data)
    f_gstar2.splt <- is.VecMatDf(f_gstar2) && NROW(as.data.frame(f_gstar2)) == NROW(data)
    if (!is.null(f_gstar2)) { # Matrices for estimators under f_gstar2
      est.communities_gstar2 <- est.communities_gstar1
      wts.communities_gstar2 <- wts.communities_gstar1
      fWi.communities_gstar2 <- fWi.communities_gstar1
      QY.communities_gstar2 <- QY.communities_gstar1
    }
    
    for (i in 1:length(communityList)) {
      if (verbose) {
        message("###########################################################################")
        message("Fitting TMLE on the " %+% i %+% "th community: " %+% communityList[i])
        message("###########################################################################")
      }
      
      ## Create an R6 object that stores and manages the subdata for each community, later passed on to estimation algorithm(s)
      subdata <- data[(data[, communityID] == communityList[i]), ]
      sub.obs.wts <- obs.wts[data[, communityID] == communityList[i]]
      if (f_gstar1.splt) { sub.f_gstar1 <- f_gstar1[data[, communityID] == communityList[i]] } else { sub.f_gstar1 <- f_gstar1 }
      if (f_gstar2.splt) { sub.f_gstar2 <- f_gstar2[data[, communityID] == communityList[i]] } else { sub.f_gstar2 <- f_gstar2 }
      inputYs <- CreateInputs(subdata[, Ynode], Qbounds, alpha, maptoYstar)
      subdata[, Ynode] <- inputYs$Ystar
      OData.ObsP0 <- DatKeepClass$new(Odata = subdata, nodes = nodes, norm.c.sVars = FALSE)
      OData.ObsP0$addYnode(YnodeVals = inputYs$Ystar)
      OData.ObsP0$addObsWeights(obs.wts = sub.obs.wts)
      nobs <- OData.ObsP0$nobs
      if (is.null(YnodeDet)) {
        determ.Q <- rep_len(FALSE, nobs)
      } else {
        determ.Q <- (data[, YnodeDet] == 1)
      }
      if (length(unique(sub.obs.wts)) > 1 && any(unlist(OData.ObsP0$type.sVar[Anodes]) != "binary")) {
        warning("sub.obs.wts are currently implemented on binary A. The results for non-binary A with weights may be unrealiable.")
      }
      
      #----------------------------------------------------------------------------------
      # Defining and estimating outcome mechanism E(Y|A, E, W)
      #----------------------------------------------------------------------------------
      if (verbose) {
        message("================================================================")
        message("fitting E(Y|A,W,E):= ", "P(" %+% nodes$Ynode %+% "=1 | " %+% paste(Q.sVars$predvars, collapse = ",") %+% ")")
        message("================================================================")
      }
      Qreg <- RegressionClass$new(outvar = Q.sVars$outvars, 
                                  predvars = Q.sVars$predvars, 
                                  subset_vars = !determ.Q, 
                                  estimator = getopt("Qestimator"))
      model.Q.init <- BinaryOutModel$new(reg = Qreg)$fit(overwrite = FALSE, data = OData.ObsP0, savespace = TRUE)
    
      #----------------------------------------------------------------------------------
      # Create an list with model estimates, data & other information that is passed on to treatment estimation procedure
      #----------------------------------------------------------------------------------
      estinfo_list <- list(
        data = subdata, 
        nodes = nodes,
        communityID = communityID,
        community.step = community.step,
        pooled.Q = pooled.Q,
        TMLE.targetStep = TMLE.targetStep,
        obs.wts = sub.obs.wts, 
        community.wts = community.wts,
        Qbounds = inputYs$Qbounds,
        lbound = lbound,
        merged.form = merged.form, 
        model.Q.init = model.Q.init,
        Q.sVars = Q.sVars,
        f.g0 = f_g0,
        h.g0.sVars = h.g0.sVars,
        h.gstar.sVars = h.gstar.sVars,
        h.g0_GenericModel = h.g0_GenericModel,
        h.gstar_GenericModel = h.gstar_GenericModel,
        savetime.fit.hbars = savetime.fit.hbars,
        n_MCsims = n_MCsims
      ) 
      estinfo_list_g1 <- append(estinfo_list, list(f.gstar = sub.f_gstar1))
      if (!is.null(f_gstar2)) { estinfo_list_g2 <- append(estinfo_list, list(f.gstar = sub.f_gstar2)) }
      
      #----------------------------------------------------------------------------------
      # Running MC evaluation for substitution TMLE, MLE and IPTW estsimators
      #----------------------------------------------------------------------------------
      # Incl. estimate treatment mechanism f(a|E, W)) and clever covariates & targeting step
      tmle_gstar1_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g1)
      if (!is.null(f_gstar2)) {
        tmle_gstar2_out <- CalcAllEstimators(OData.ObsP0 = OData.ObsP0, est_params_list = estinfo_list_g2)
      } else {
        tmle_gstar2_out <- NULL
      }
      
      #----------------------------------------------------------------------------------
      # Store substitution estsimators for each corresponding community
      #----------------------------------------------------------------------------------
      est.communities_gstar1[i, ] <- tmle_gstar1_out$ests_mat[, 1]
      wts.communities_gstar1 <- rbind(wts.communities_gstar1, tmle_gstar1_out$wts_mat)
      fWi.communities_gstar1 <- rbind(fWi.communities_gstar1, tmle_gstar1_out$fWi_mat)
      QY.communities_gstar1 <- rbind(QY.communities_gstar1, tmle_gstar1_out$QY_mat)
      obs.wts.communities <- c(obs.wts.communities, tmle_gstar1_out$obs.wts)
      if (!is.null(f_gstar2)) {
        est.communities_gstar2[i, ] <- tmle_gstar2_out$ests_mat[, 1]
        wts.communities_gstar2 <- rbind(wts.communities_gstar2, tmle_gstar2_out$wts_mat)
        fWi.communities_gstar2 <- rbind(fWi.communities_gstar2, tmle_gstar2_out$fWi_mat)
        QY.communities_gstar2 <- rbind(QY.communities_gstar2, tmle_gstar2_out$QY_mat)
      }
    }
    
    # Reconstruct the results of substitution estsimators based on the combined matrices
    community.wts.pair <- community.wts[match(communityList, community.wts[, "id"]), "weights"]
    est_mat_gstar1 <- matrix(0L, nrow = 3, ncol = 1)
    est_mat_gstar1[, 1] <- apply(est.communities_gstar1, 2, weighted.mean, w = community.wts.pair)
    rownames(est_mat_gstar1) <- c("TMLE", "IPTW", "MLE"); colnames(est_mat_gstar1) <- "estimate"
    tmle_gstar1.communities <- list(ests_mat = est_mat_gstar1, wts_mat = wts.communities_gstar1, fWi_mat = fWi.communities_gstar1, 
                                    QY_mat = QY.communities_gstar1, obs.wts = obs.wts.communities)
    if (!is.null(f_gstar2)) {
      est_mat_gstar2 <- matrix(0L, nrow = 3, ncol = 1)
      est_mat_gstar2[, 1] <- apply(est.communities_gstar2, 2, weighted.mean, w = community.wts.pair)
      rownames(est_mat_gstar2) <- c("TMLE", "IPTW", "MLE"); colnames(est_mat_gstar2) <- "estimate"
      tmle_gstar2.communities <- list(ests_mat = est_mat_gstar2, wts_mat = wts.communities_gstar2, fWi_mat = fWi.communities_gstar2, 
                                      QY_mat = QY.communities_gstar2, obs.wts = obs.wts.communities)
    }
    
    #----------------------------------------------------------------------------------
    # Create output list (estimates, as. variances, CIs)
    #----------------------------------------------------------------------------------
    inputYs <- CreateInputs(data[, Ynode], Qbounds, alpha, maptoYstar)
    data[, Ynode] <- inputYs$Ystar
    OData.ObsP0 <- DatKeepClass$new(Odata = data, nodes = nodes, norm.c.sVars = FALSE)
    OData.ObsP0$addYnode(YnodeVals = inputYs$Ystar)
    OData.ObsP0$addObsWeights(obs.wts = obs.wts)
    tmle_gstar1.communities$OData.ObsP0 <- OData.ObsP0
    if (!is.null(f_gstar2)) { tmle_gstar2.communities$OData.ObsP0 <- OData.ObsP0 }
    
    # **** Double check IC-based variance calculation with Prof Mark ****
    EY_gstar1 <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, tmle_g_out = tmle_gstar1.communities)
    EY_gstar2 <- NULL
    ATE <- NULL	
    if (!is.null(f_gstar2)) {
      EY_gstar2 <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, tmle_g_out = tmle_gstar2.communities)
      ATE <- calcParameters(inputYs = inputYs, alpha = CI_alpha, est_params_list = estinfo_list, 
                            tmle_g_out = tmle_gstar1.communities, tmle_g2_out = tmle_gstar2.communities)
    }
  }  
  tmleCommunity.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)
  class(tmleCommunity.res) <- c(class(tmleCommunity.res), "tmleCommunity")
  return(tmleCommunity.res)
}
