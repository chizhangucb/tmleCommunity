# @title tmleCommunity-package
# @docType package
# @author Chi Zhang, Oleg Sofrygin, Mark J. van der Laan

#' @useDynLib tmleCommunity
#' @import R6
#' @importFrom grDevices nclass.FD nclass.Sturges nclass.scott
#' @importFrom graphics axis barplot hist par text
#' @importFrom methods is
#' @importFrom stats approx quasibinomial binomial coef glm.control glm.fit plogis predict qlogis qnorm quantile rnorm runif terms var
#' @importFrom utils data head str
#' @importFrom Hmisc wtd.var

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
  rownames(teststat_mat) <- c("tmle", "iptw", "gcomp")
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
  working.model <- est_params_list$working.model
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
  if ((community.step == "individual_level" && working.model == TRUE) || community.step == "perCommunity") { 
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
  working.model <- est_params_list$working.model
  TMLE.targetStep <- est_params_list$TMLE.targetStep
  nodes <- OData.ObsP0$nodes
  Y <- OData.ObsP0$noNA.Ynodevals # actual observed & transformed Y's
  determ.Q <- OData.ObsP0$det.Y
  model.Q.init <- est_params_list$model.Q.init
  
  # getting all node vals, including deterministic Y and predictions P(Y=1) for non-DET Y under orignal A
  QY.init <- OData.ObsP0$noNA.Ynodevals 
  QY.init[!OData.ObsP0$det.Y] <- model.Q.init$predict(newdata = OData.ObsP0)$getprobA1[!OData.ObsP0$det.Y]
  QY.init <- bound(QY.init, est_params_list$Qbounds)
  off <- qlogis(QY.init)  # offset
  
  if (community.step == "individual_level" && working.model == FALSE) { # if we do NOT believe our working model (i.e. estimate under the lareg model)
    # if (!is.null(communityID)) {}  # Don't need it here since tmleCommunity takes care of it
    # Since individual-level TMLE with no working.model requires 'communityID' to aggregate data to the cluster-level in the estimation of 
    # trt mechanism. Lack of 'communityID' pool data over all communities & treat it as non-hierarchical data when fitting clever covariates
    
    # aggregate initial outcome predictions to the cluster-level & Recalculate offset based on aggregated initiate predictions
    QY.init <- aggregate(x = QY.init, by=list(id = data[, communityID]), mean)[, 2]
    off <- qlogis(QY.init)  
    # aggregate original dataset to the cluster-level & redefine OData.ObsP0
    Y <- aggregate(x = Y, by=list(id = data[, communityID]), mean)[, 2] 
    determ.Q <- rep_len(FALSE, length(Y))  # For aggregated data, YnodeDet is currently unavailable, treat all Y^c as nondeterministic
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
  obs.wts <- est_params_list$obs.wts <- OData.gstar$get.obsweights(TRUE)
  
  #************************************************
  # IPTW_h estimator:
  #************************************************  
  IPTW <- Y
  IPTW[!determ.Q] <- Y[!determ.Q] * h_wts[!determ.Q]
  # IPTW_unwt <- mean(IPTW)
  if (community.step == "individual_level" && working.model == TRUE) {
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
#' Estimate Marginal Treatment Effects For Arbitrary (Stochastic) Interventions
#'
#' Estimate the marginal treatment effect among i.i.d units using \strong{TMLE} (targeted maximum likelihood estimation). It also provide
#' \strong{IPTW} (the inverse-probability-of-treatment or Horvitz-Thompson) and \strong{GCOMP} (parametric G-computation formula).
#' @param data \code{data.frame} with named columns, containing \code{WEnodes}, \code{Anode}, \code{Ynode} and possibly \code{communityIndex}.
#' @param Ynode Column names or indices in \code{data} of outcome variable name. Outcome can be either binary or continuous. 
#'   This can instead be specified on the left-side of the regression formula in argument \code{Qform}.
#' @param Anodes Column names or indices in \code{data} of exposure (treatment) variables; exposures can be either binary, categorical or continuous.
#' @param WEndoes Column names or indices in \code{data} of individual-level (and possibly community-level) baseline covariates.
#'   Factors are not currently allowed.
#' @param communityID Optional column name or index in \code{data} of community identifier variable. If known, either stratify on community level
#'   when estimating outcome and treatment mechanisms, or perform panel transformation on data before estiamtion, depending on \code{community.step}.
#' @param YnodeDet Optional column name or index in \code{data} of deterministic values of outcome Ynode, coded as (TRUE/FALSE) or (1/0). If TRUE/1, 
#'  value of Ynode is given deterministically / constant. 
#' @param community.step Methods to deal with community-level data, one of "NoCommunity" (Default), "community_level" and "individual_level". 
#'  If communityID = NULL, then automatically pool over all communities.
#' @param working.model Logical
#' @param community.wts Optional matrix of community-level observation weights (where dim = the number of communities by 2). The first   
#'  column contains the communities' names (ie., \code{data[, communityID]}) and the second column contains the corresponding weights.   
#'  If "equal.community", assumed to be all 1. Currently only support a numeric vector, "equal.community" (Default) and "size.community" 
#'  (i.e., by setting community.wts = "size.community", treat the number of individuals within each community as its weight, respectively).
#' @param f_gstar1 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  If a matrix/ data frame, its number of rows must be either nrow(data) or 1 (constant exposure assigned to all observations), and its number of 
#'  columns must be length(Anodes). If a vector, it must be of length nrow(data) or 1. If a function, it must return a data frame of counterfactual
#'  exposures sampled based on Anodes, WEnodes (and possibly communityIndex) passed as a named argument "data". Thus, the function must 
#'  include "data" as one of its argument names. The interventions defined by f_gstar1 can be static, dynamic or stochastic. See Exmaples below.
#' @param f_gstar2 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  It has the same components and requirements as f_gstar1
#' @param Qform Character vector of regression formula for Ynode. If not specified, the outcome variable is regressed on all covariates included in 
#'  Anodes and WEnodes.
#' @param Qbounds Upper and lower bounds on Y and predicted values for initial Q. Defaults to the range of Y, widened by 10\% of the min and max values.
#' @param alpha Used to keep predicted values for initial Q bounded away from (0,1) for logistic fluctuation.
#' @param fluctuation Default to "logistic", it could also be "linear" (for targeting step).
#' @param f_g0 Optional function used to specify model knowledge about value of Anodes. It estimates \code{P(A | W, E)} under \code{g0} by 
#'  sampling a large vector/ data frame of Anode (of length or number of rows \code{nrow(data)*n_MCsims}) from \code{f_g0}
#' @param hform.g0 Character vector of regression formula for estimating the conditional density of P(A | W, E) under the observed treatment mechanism
#'  g0. If not specified, its form will be Anodes ~ WEnodes. If there are more than one expsosure, it fits a joint probability.
#' @param hform.gstar Character vector of regression formula for estimating the conditional density P(A | W, E) under interventions f_gstar1 or f_gstar2. 
#'  If not specified, it follows the same rule used in hform.g0. 
#' @param lbound Value between (0,1) for truncation of predicted P(A | W, E). Default to 0.005
#' @param obs.wts Optional vector of individual-level observation (sampling) weights (of length \code{nrow(data)}). If NULL, assumed to be all 1. 
#' @param h.g0_GenericModel Previously fitted models for P(A | W, E) under g0, one of returns of tmleCommunity.function. If known, predictions
#'  for P(A=a | W=w, E=e) under g0 are based on the fitted models in \code{h.g0_GenericModel}.
#' @param h.gstar_GenericModel Previously fitted models for P(A^* | W, E) under gstar, one of returns of tmleCommunity.function. If known,  
#'  Predictions for P(A=a | W=w, E=e) under gstar are based on the fitted models in \code{h.gstar_GenericModel}.
#' @param savetime.fit.hbars Logical for saving time when fitting fit.hbars. If true, 
#' @param TMLE.targetStep TMLE targeting step method, either "tmle.intercept" (Default) or "tmle.covariate".
#' @param n_MCsims Number of simulations for Monte-Carlo analysis. Each simulation generates new exposures under f_gstar1 or f_gstar2, with a sample 
#'  size of nrow(data). Then these expsosures are used when fitting the conditional densities P(A | W, E).
#' @param rndseed Random seed for controlling sampling A under f_gdelta1 or f_gdelta2 (reproducibility)
#' @param verbose Flag. If TRUE, print status messages. Default to TRUE.
#'
#' @section IPTW estimator:
#' **********************************************************************
#'
#' @section GCOMP estimator:
#' **********************************************************************
#'
#' @section TMLE estimator:
#' **********************************************************************
#'
#' @section Modeling \code{P(A|W)} for covariates \code{(A,W)}:
#' **********************************************************************
#' Non-parametric
#'  estimation of the common \strong{unit-level} multivariate joint conditional probability model \code{P_g0(A|W)},
#'  for unit-level summary measures \code{(sA,sW)} generated from the observed exposures and baseline covariates
#'  \eqn{(A,W)=(A_i,W_i : i=1,...,N)} (their joint density given by \eqn{g_0(A|W)Q(W)}), is performed by first
#'  constructing the dataset of N summary measures, \eqn{(sA_i,sW_i : i=1,...,N)}, and then fitting the usual i.i.d. MLE
#'  for the common density \code{P_g0(A|W)} based on the pooled N sample of these summary measures.
#'
#'  Note that \code{A} can be multivariate and any of its components \code{A[j]} can be either binary, categorical
#'  or continuous. The joint probability model for \code{P(A|W)} = \code{P(A[1],...,A[k]|W)} can be factorized as
#'  \code{P(A[1]|W)} * \code{P(A[2]|W, A[1])} * ... * \code{P(A[k]|W, A[1],...,A[k-1])},
#'  where each of these conditional probability models is fit separately, depending on the type of the outcome variable \code{A[j]}.
#'
#'  If \code{A[j]} is binary, the conditional probability \code{P(A[j]|W,A[1],...,A[j-1])} is evaluated via logistic regression model.
#'  When \code{sA[j]} is continuous (or categorical), its range will be fist partitioned into \code{K} bins and the corresponding \code{K}
#'  bin indicators (\code{B_1,...,B_K}), where each bin indicator \code{B_j} is then used as an outcome in a
#'  separate logistic regression model with predictors given by \code{W, A[1],..., A[k-1]}.
#'  Thus, the joint probability \code{P(A|W)} is defined by such a tree of binary logistic regressions.
#'
#' For simplicity, we now suppose \code{sA} is continuous and univariate and we describe here an algorithm for fitting \eqn{P_{g_0}(A|W)} 
#'  (the algorithm for fitting \eqn{P_{g^*}(A^*|W^*)} is equivalent, except that exposure \code{A} is replaced with exposure \code{A^*}
#'  generated under \code{f_gstar1} or \code{f_gstar2} and the predictors \code{W} from the regression formula
#'  \code{hform.g0} are replaced with predictors \code{W^*} specified by the regression formula \code{hform.gstar}).
#'
#' \enumerate{
#' \item Generate a dataset of N observed continuous summary measures (\code{A_i}:i=1,...,N) from observed ((\code{A_i},\code{W_i}):i=1,...,N).
#'
#' \item Divide the range of \code{sA} values into intervals S=(i_1,...,i_M,i_{M+1}) so that any observed data point
#'    \code{sa_i} belongs to one interval in S, namely, for each possible value sa of \code{A} there is k\\in{1,...,M}, such that, 
#'    i_k < \code{a} <= i_{k+1}. Let the mapping B(sa)\\in{1,...,M} denote a unique interval in S for a, such that, i_{B(a)} < a <= i_{B(a)+1}.
#'    Let bw_{B(a)}:=i_{B(a)+1}-i_{B(a)} be the length of the interval (bandwidth) (i_{B(a)},i_{B(a)+1}).
#'    Also define the binary indicators b_1,...,b_M, where b_j:=I(B(a)=j), for all j <= B(a) and b_j:=NA for all j>B(a). That is we set b_j to 
#'    missing ones the indicator I(B(a)=j) jumps from 0 to 1. Now let \code{A} denote the random variable for the observed exposure for one unit
#'    and denote by (B_1,...,B_M) the corresponding random indicators for \code{A} defined as B_j := I(B(\code{A}) = j) for all j <= B(\code{A}) 
#'    and B_j:=NA for all j>B(\code{A}).
#'
#' \item For each j=1,...,M, fit the logistic regression model for the conditional probability P(B_j = 1 | B_{j-1}=0, W), i.e.,
#'    at each j this is defined as the conditional probability of B_j jumping from 0 to 1 at bin j, given that B_{j-1}=0 and
#'    each of these logistic regression models is fit only among the observations that are still at risk of having B_j=1 with B_{j-1}=0.
#'
#' \item Normalize the above conditional probability of B_j jumping from 0 to 1 by its corresponding interval length (bandwidth) bw_j to
#'    obtain the discrete conditional hazards h_j(W):=P(B_j = 1 | (B_{j-1}=0, W) / bw_j, for each j.
#'    For the summary measure \code{A}, the above conditional hazard h_j(sW) is equal to P(\code{A} \\in (i_j,i_{j+1}) | \code{A}>=i_j, sW),
#'    i.e., this is the probability that \code{A} falls in the interval (i_j,i_{j+1}), conditional on sW and conditional on the fact that
#'    \code{A} does not belong to any intervals before j.
#'
#' \item  Finally, for any given data-point \code{(a,w)}, evaluate the discretized conditional density for P(\code{A}=a|W=w) by first
#'    evaluating the interval number k=B(a)\\in{1,...,M} for \code{a} and then computing \\prod{j=1,...,k-1}{1-h_j(W))*h_k(W)}
#'    which is equivalent to the joint conditional probability that \code{a} belongs to the interval (i_k,i_{k+1}) and does not belong
#'    to any of the intervals 1 to k-1, conditional on sW.
#'  }
#'
#' The evaluation above utilizes a discretization of the fact that any continuous density f of random variable X can be written as f_X(x)=S_X(x)*h_X(x),
#'  for a continuous density f of X where S_X(x):=P(X>x) is the survival function for X, h_X=P(X>x|X>=x) is the hazard function for X; as well as the fact that
#'  the discretized survival function S_X(x) can be written as a of the hazards for s<x: S_X(x)=\\prod{s<x}h_X(x).
#'
#' @section Three methods for defining bin (interval) cuttoffs for a continuous one-dimenstional summary measure \code{A[j]}:
#' **********************************************************************
#'
#' There are 3 alternative methods to defining the bin cutoffs S=(i_1,...,i_M,i_{M+1}) for a continuous summary measure
#'  \code{A}. The choice of which method is used along with other discretization parameters (e.g., total number of
#'  bins) is controlled via the tmlenet_options() function. See \code{?tmlenet_options} argument \code{bin.method} for
#'  additional details.
#'
#' Approach 1 (\code{equal.len}): equal length, default.
#'
#' *********************
#'
#' The bins are defined by splitting the range of observed \code{A} (sa_1,...,sa_n) into equal length intervals.
#'  This is the dafault discretization method, set by passing an argument \code{bin.method="equal.len"} to
#'  \code{tmlenet_options} function prior to calling \code{tmleCommunity()}. The intervals will be defined by splitting the
#'  range of (sa_1,...,sa_N) into \code{nbins} number of equal length intervals, where \code{nbins} is another argument
#'  of \code{tmleCom_Options()} function. When \code{nbins=NA} (the default setting) the actual value of \code{nbins}
#'  is computed at run time by taking the integer value (floor) of \code{n/maxNperBin},
#'  for \code{n} - the total observed sample size and \code{maxNperBin=1000} - another argument of
#'  \code{tmleCom_Options()} with the default value 1,000.
#'
#' Approach 2 (\code{equal.mass}): data-adaptive equal mass intervals.
#'
#' *********************
#'
#' The intervals are defined by splitting the range of \code{A} into non-equal length data-adaptive intervals that
#'  ensures that each interval contains around
#'  \code{maxNperBin} observations from (sa_j:j=1,...,N).
#'  This interval definition approach can be selected by passing an argument \code{bin.method="equal.mass"} to
#'  \code{tmleCom_Options()} prior to calling \code{tmleCommunity()}.
#'  The method ensures that an approximately equal number of observations will belong to each interval, where that number
#'  of observations for each interval
#'  is controlled by setting \code{maxNperBin}. The default setting is \code{maxNperBin=1000} observations per interval.
#'
#' Approach 3 (\code{dhist}): combination of 1 & 2.
#'
#' *********************
#'
#' The data-adaptive approach dhist is a mix of Approaches 1 & 2. See Denby and Mallows "Variations on the Histogram"
#'  (2009)). This interval definition method is selected by passing an argument \code{bin.method="dhist"} to
#'  \code{tmleCom_Options()}  prior to calling \code{tmleCommunity()}.
#'
#' @return A named list with 3 items containing the estimation results for:
#'  \itemize{
#'  \item \code{EY_gstar1} - estimates of the mean counterfactual outcome under (stochastic) intervention function \code{f_gstar1} \eqn{(E_{g^*_1}[Y])}.
#'  \item \code{EY_gstar2} - estimates of the mean counterfactual outcome under (stochastic) intervention function \code{f_gstar2} \eqn{(E_{g^*_2}[Y])}, 
#'    or \code{NULL} if \code{f_gstar2} not specified.
#'  \item \code{ATE} - additive treatment effect (\eqn{E_{g^*_1}[Y]} - \eqn{E_{g^*_2}[Y]}) under interventions \code{f_gstar1}
#'    vs. in \code{f_gstar2}, or \code{NULL} if \code{f_gstar2} not specified.
#' }
#' Each list item above is itself a list containing the items:
#'  \itemize{
#'  \item \code{estimates} - various estimates of the target parameter (network population counterfactual mean under (stochastic) intervention).
#'  \item \code{vars} - the asymptotic variance estimates, for \strong{TMLE}, \strong{IPTW} and \strong{GCOMP}. Notice, inference for gcomp is 
#'    not accurate! It is based on TMLE influence curves.
#'  \item \code{CIs} - CI estimates at \code{alpha} level, for \strong{TMLE}, \strong{IPTW} and \strong{GCOMP}.
#'  \item \code{h.g0_GenericModel} - The model fits for P(\code{A}|\code{W, E}) under observed exposure mechanism
#'    \code{g0}. This is an object of \code{h.g0_GenericModel} \pkg{R6} class.
#'  \item \code{h.gstar_GenericModel} - The model fits for P(\code{A}|\code{W, E}) under intervention \code{f_gstar1}
#'    or \code{f_gstar2}. This is an object of \code{GenericModel} \pkg{R6} class.
#' }
#' Currently implemented estimators are:
#'  \itemize{
#'  \item \code{tmle} - Either weighted regression intercept-based TMLE (\code{tmle.intercept} - the default) with weights defined by the IPTW weights
#'    \code{h_gstar/h_gN} or covariate-based unweighted TMLE (\code{tmle.covariate}) that uses the IPTW weights as a covariate \code{h_gstar/h_gN}.  
#'  \item \code{iptw} - Efficient IPTW based on weights h_gstar/h_gN.
#'  \item \code{gcomp} - Parametric G-computation formula substitution estimator.
#' }
#' @example tests/examples/3_tmleCommunity_examples.R
#' @export
tmleCommunity <- function(data, Ynode, Anodes, WEnodes, YnodeDet = NULL, obs.wts = c("equal.within.pop", "equal.within.community"), 
                          community.wts = c("equal.community", "size.community"), communityID = NULL, working.model = FALSE, 
                          community.step = c("NoCommunity", "community_level", "individual_level", "perCommunity"), 
                          f_g0 = NULL, f_gstar1, f_gstar2 = NULL, Qform = NULL, Qbounds = NULL, alpha = 0.995,                                                      
                          fluctuation = "logistic", hform.g0 = NULL, hform.gstar = NULL, lbound = 0.005, 
                          h.g0_GenericModel = NULL, h.gstar_GenericModel = NULL, savetime.fit.hbars = TRUE, 
                          TMLE.targetStep = c("tmle.intercept", "tmle.covariate"),
                          n_MCsims = 1, CI_alpha = 0.05, rndseed = NULL, verbose = TRUE) {
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  gvars$verbose <- verbose
  if (verbose) message("Running tmleCommunity with the following settings from tmleCom_Options(): "); str(gvars$opts)
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (is.null(savetime.fit.hbars)) savetime.fit.hbars <- getopt("savetime.fit.hbars")
  if (obs.wts == "equal.within.pop") { # weigh individuals in the entire dataset equally so big community gets bigger total weight
    obs.wts <- rep(1, NROW(data))
  } else if (obs.wts == "equal.within.community") { # weigh individuals in each community equally and weigh communities equally
    obs.wts <- rep(as.vector(1/table(data[, communityID])), as.vector(table(data[, communityID])))
  }
  if (is.character(community.step) && (community.step != "NoCommunity") && !is.data.frame(community.wts)) {
    community.wts.mat <- as.data.frame(matrix(0L, nrow = length(unique(data[, communityID])), ncol = 2))
    colnames(community.wts.mat) <- c("id", "weights")
    community.wts.mat[, 1] <- names(table(data[, communityID]))
    if (community.wts == "equal.community") { # weigh each community equally
      community.wts.mat[, 2]  <- rep(1, length(unique(data[, communityID]))) 
    } else if (community.wts == "size.community") { # weigh each community by its number of observations - The larger community has larger weight
      community.wts.mat[, 2] <- as.vector(table(data[, communityID]))
    } 
    community.wts <- community.wts.mat; community.wts.mat <- NULL
  }
  if (!is.null(Qform) && !is.null(Ynode)) {
    Qform <- paste(Ynode, substring(Qform, first = as.numeric(gregexpr("~", Qform))))
    message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+% Ynode)
    message("Thus the Qform becomes " %+% Qform)
  }
  if (!is.null(Qform) && is.null(Ynode)) {
    Ynode <- LhsVars(Qform)[1]
    message("Setting the Ynode to: " %+% Ynode)
  }
  community.step <- community.step[1]
  TMLE.targetStep <- TMLE.targetStep[1]
  
  ## Check if any unexpected inputs
  if (!(TMLE.targetStep %in% c("tmle.intercept", "tmle.covariate"))) 
    stop("TMLE.targetStep argument must be either 'tmle.intercept' or 'tmle.covariate'")
  if (!(community.step %in% c("NoCommunity", "community_level", "individual_level", "perCommunity"))) 
    stop("community.step argument must be one of 'NoCommunity', 'community_level', 'individual_level' and 'perCommunity'")
  if ((community.step %in% c("community_level", "individual_level", "perCommunity")) & is.null(communityID)) {
    messageMSg <- c("'communityID' is required when using 'community_level', 'individual_level' and 'perCommunity. '",
                    "Lack of 'communityID' forces the algorithm to automatically pool data over all communities ",
                    "and treat it as non-hierarchical dataset. Details in package documentation.\n",
                    "In other words, we simply treat community.step = 'NoCommunity' and working.model = FALSE")
    message(messageMSg[1] %+% messageMSg[2] %+% messageMSg[3] %+% messageMSg[4])
    community.step <- "NoCommunity"; working.model <- FALSE
  } 
  if (!(community.wts %in% c("equal.community", "size.community")) && !is.data.frame(community.wts)) {
    stop("Currently only numeric values, 'equal.community' and 'size.community' are supported for community.wts")
  }
  nodes <- list(Ynode = Ynode, Anodes = Anodes, WEnodes = WEnodes, communityID = communityID)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, Qform, hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts, community.wts)) stop()
  colnames(community.wts) <- c("id", "weights")  # ensure that the column names are defined as "id" and "weights"
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
  if (community.step == "community_level") { # if running entire TMLE algorithm at cluster-level, aggregate data now
    # if (!is.null(communityID)) {
    data <- aggregate(x = data, by=list(newid = data[, communityID]), mean) # [, 2 : (ncol(data)+1)] # Don't keep the extra ID column
    colname.allNA <- colnames(data)[colSums(is.na(data)) == NROW(data)]  # columns with all NAs after aggregation, due to non-numeric values
    if (length(colname.allNA) != 0) {
      data <- data[, -which(colnames(data) %in% colname.allNA)] # Remove columns that contain only NAs inside
      names(data)[which(colnames(data) == "newid")] <- communityID  # change 'newid' back to the communityID name
      warning(paste(colname.allNA, collapse = ', ') %+% " is(are) removed from the aggregated data due to all NAs in the column(s).")
      warning("Suggestion: convert the non-numeric values to numeric, e.g., create dummy variables for each category/ string/ factor.")
    }
    obs.wts <- community.wts[match(data[, communityID], community.wts[, "id"]), "weights"]  # ensure that weights match with their corresponding communities
    # } else {
    #   warningMesg <- c("Since community-level TMLE requires communityID to aggregate to the cluster-level. Lack of 'communityID' forces the ",
    #                    "algorithm to automatically pool data over all communities and treat it as non-hierarchical dataset")
    #   warning(warningMesg[1] %+% warningMesg[2])
    # }
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
      warning("obs.wts are currently implemented on binary A. The results for non-binary A with weights may be unrealiable.")
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
      working.model = working.model,
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
    message("######################################################################################")
    message("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.")
    message("######################################################################################")
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
    if (!is.null(f_gstar2)) { # Matrices for estimators under f_gstar2
      est.communities_gstar2 <- est.communities_gstar1
      wts.communities_gstar2 <- wts.communities_gstar1
      fWi.communities_gstar2 <- fWi.communities_gstar1
      QY.communities_gstar2 <- QY.communities_gstar1
    }
    
    for (i in 1:length(communityList)) {
      message("###########################################################################")
      message("Fitting TMLE on the " %+% i %+% "th community: " %+% communityList[i])
      message("###########################################################################")
      
      ## Create an R6 object that stores and manages the subdata for each community, later passed on to estimation algorithm(s)
      subdata <- data[(data[, communityID] == communityList[i]), ]
      sub.obs.wts <- obs.wts[data[, communityID] == communityList[i]]
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
        working.model = working.model,
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
      estinfo_list_g1 <- append(estinfo_list, list(f.gstar = f_gstar1))
      if (!is.null(f_gstar2)) { estinfo_list_g2 <- append(estinfo_list, list(f.gstar = f_gstar2)) }
      
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
