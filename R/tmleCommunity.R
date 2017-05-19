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
calcParameters <- function(OData.ObsP0, inputYs, alpha = 0.05, tmle_g_out, tmle_g2_out = NULL) {
  nobs <- OData.ObsP0$nobs
  ests_mat <- tmle_g_out$ests_mat
  QY_mat <- tmle_g_out$QY_mat
  fWi_mat <- tmle_g_out$fWi_mat
  wts_mat <- tmle_g_out$wts_mat
  obs.wts <- tmle_g_out$obs.wts
  # ests_unwt_mat <- tmle_g_out$ests_unwt_mat
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
  var_mat.res <- get_est_sigmas(estnames = c("tmle", "iptw", "gcomp"), obsYvals = OData.ObsP0$noNA.Ynodevals, 
                                ests_mat = ests_mat, QY_mat = QY_mat, wts_mat = wts_mat, fWi_mat = fWi_mat, obs.wts = obs.wts)
  as.var_mat <- var_mat.res$as.var_mat
  if (maptoYstar) {
    as.var_mat <- as.var_mat * (diff(ab) ^ 2)
    if (is.null(tmle_g2_out)) {
      ests_mat <- ests_mat * diff(ab) + ab[1]
    } else { 
      ests_mat <- ests_mat * diff(ab) 
    }
  }
  
  get_CI <- function(xrow, n) {
    f_est_CI <- function(n, psi, sigma2_N) { # get CI
      z_alpha <- qnorm(1-alpha/2)
      CI_est <- c(psi - z_alpha * sqrt(sigma2_N / n), psi + z_alpha * sqrt(sigma2_N / n))
      return(CI_est)
    }
    psi <- xrow["estimate"];
    sigma2_N <- xrow["Var"];
    return(f_est_CI(n = n, psi = psi, sigma2_N = sigma2_N))
  }
  
  CIs_mat <- t(apply(cbind(ests_mat, as.var_mat), 1, get_CI, n = nobs))
  colnames(CIs_mat) <- c("LBCI_" %+% as.character(alpha/2), "UBCI_" %+% as.character(1-alpha/2))
  
  # ----------------------------------------------------
  # RENAME ESTIMATORS FOR THE FINAL OUTPUT:
  # ----------------------------------------------------
  rownames(ests_mat) <- c("tmle", "iptw", "gcomp")
  rownames(as.var_mat) <- c("tmle", "iptw", "gcomp")
  rownames(CIs_mat) <- c("tmle", "iptw", "gcomp")
  
  EY_g.star <- list(estimates = ests_mat, 
                    # est_unwt = ests_unwt_mat,
                    vars = (as.var_mat / nobs), 
                    CIs = CIs_mat, 
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
get_est_sigmas <- function(estnames, obsYvals, ests_mat, QY_mat, wts_mat, fWi_mat, obs.wts) {
  fWi <- fWi_mat[, "fWi_Qinit"]
  QY.init <- QY_mat[, "QY.init"] 
  h_wts <- wts_mat[, "h_wts"]
  
  # TMLE inference based on the iid IC: (** Use QY.init not QY.star)
  iidIC_tmle <- h_wts * (obsYvals - QY.init) + fWi - ests_mat[rownames(ests_mat) %in% "TMLE",]
  var_iid.tmle <- wtd.var(iidIC_tmle, weights = obs.wts, normwt = T)
  # var_iid.tmle <- mean((iidIC_tmle)^2)  # assume mean(iidIC_tmle) = 0
  
  # MLE inference based on the iid IC: (** Use QY.init not QY.star) *** NOT ACCURATE
  iidIC_mle <- h_wts * (obsYvals - QY.init) + fWi - ests_mat[rownames(ests_mat) %in% "MLE",]
  var_iid.mle <- wtd.var(iidIC_mle, weights = obs.wts, normwt = T)
  # var_iid.mle <- mean((iidIC_mle)^2)  # assume mean(iidIC_mle) = 0
  
  # IPTW h (based on the mixture density clever covariate (h)):
  iidIC_iptw <- h_wts * (obsYvals) - ests_mat[rownames(ests_mat) %in% "IPTW",]
  var_iid.iptw <- wtd.var(iidIC_iptw, weights = obs.wts, normwt = T)
  # var_iid.iptw <- mean((iidIC_iptw)^2)  # assume mean(iidIC_iptw) = 0
  
  as.var_mat <- matrix(0, nrow = 3, ncol = 1)
  as.var_mat[, 1] <- c(var_iid.tmle, var_iid.iptw, var_iid.mle)
  rownames(as.var_mat) <- estnames; colnames(as.var_mat) <- "Var"
  return(list(as.var_mat = as.var_mat, IC = list(IC.tmle = iidIC_tmle, IC.iptw = iidIC_iptw, IC.gcomp = iidIC_mle)))
}


#--------------------------------------- CalcAllEstimators ---------------------------------
# Purpose: Estimate h_bar under g0 and g* given observed data and vector of c^Y's data is an DatKeepClass object
#-------------------------------------------------------------------------------------------
CalcAllEstimators <- function(OData.ObsP0, est_params_list) {
  TMLE.targetStep <- est_params_list$TMLE.targetStep
  nodes <- OData.ObsP0$nodes
  Y <- OData.ObsP0$noNA.Ynodevals # actual observed & transformed Y's
  determ.Q <- OData.ObsP0$det.Y
  model.Q.init <- est_params_list$model.Q.init
  
  QY.init <- OData.ObsP0$noNA.Ynodevals # getting all node vals, inc. deterministic  
  QY.init[!OData.ObsP0$det.Y] <- model.Q.init$predict(newdata = OData.ObsP0)$getprobA1[!OData.ObsP0$det.Y] # predictions P(Y=1) for non-DET Y
  off <- qlogis(QY.init)  # offset
  
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
  IPTW <- weighted.mean(IPTW, w = obs.wts)
  
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
  MCS_res_wt <- MCS_res$wtmean_psis_all
  # ests_unwt <- MCS_res$unwt.mean_psis_all
  
  # ests_unwt <- c(TMLE = ests_unwt[["TMLE"]], IPTW = IPTW_unwt, MLE = ests_unwt[["MLE"]])
  # ests_unwt_mat <- matrix(0L, nrow = length(ests_unwt), ncol = 1)
  # ests_unwt_mat[, 1] <- ests_unwt
  # rownames(ests_unwt_mat) <- names(ests_unwt); colnames(ests_unwt_mat) <- "estimate"

  ests <- c(TMLE = MCS_res_wt[["TMLE"]], IPTW = IPTW, MLE = MCS_res_wt[["MLE"]])
  ests_mat <- matrix(0L, nrow = length(ests), ncol = 1)
  ests_mat[, 1] <- ests
  rownames(ests_mat) <- names(ests); colnames(ests_mat) <- "estimate"
  
  wts_mat <- matrix(0L, nrow = OData.ObsP0$nobs, ncol = 1)
  colnames(wts_mat) <- c("h_wts")
  wts_mat[, "h_wts"] <- h_wts
  
  fWi_mat <- matrix(0L, nrow = OData.ObsP0$nobs, ncol = 1)
  colnames(fWi_mat) <- c("fWi_Qinit")
  fWi_mat[,"fWi_Qinit"] <- MCS_res_wt[agrep("fWi_init_", names(MCS_res_wt), max.distance=list(insertions=0, substitutions=0))]
  
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
              # ests_unwt_mat = ests_unwt_mat, 
              obs.wts = obs.wts, 
              h.g0_GenericModel = model.h.fit$genericmodels.g0,
              h.gstar_GenericModel = model.h.fit$genericmodels.gstar))
}


#------------------------------------
#' Estimate Marginal Treatment Effects For Arbitrary (Stochastic) Interventions
#'
#' Estimate the marginal treatment effect among i.i.d units using \strong{TMLE} (targeted maximum likelihood estimation). It also provide
#' \strong{IPTW} (the inverse-probability-of-treatment or Horvitz-Thompson) and \strong{GCOMP} (parametric G-computation formula).
#' @param data \code{data.frame} with named columns, containing \code{Wnodes}, \code{Anode}, \code{Ynode} 
#'  and possibly \code{Enodes} and \code{communityIndex}.
#' @param Ynode Column names or indices in \code{data} of outcome variable name. Outcome can be either binary or continuous. 
#'   This can instead be specified on the left-side of the regression formula in argument \code{Qform}.
#' @param Anodes Column names or indices in \code{data} of exposure (treatment) variables; exposures can be either binary, categorical or continuous.
#' @param Wndoes Column names or indices in \code{data} of individual-level baseline covariates. Factors are not currently allowed.
#' @param Endoes Optional column names or indices in \code{data} of community-level baseline covariates. Factors are not currently allowed.
#' @param StratifyInd Optional column name or index in \code{data} of community identifier variable. If known, stratify on community level when 
#'  estimating outcome and treatment mechanisms. If NULL, pool over all communities. 
#' @param YnodeDet Optional column name or index in \code{data} of deterministic values of outcome Ynode, coded as (TRUE/FALSE) or (1/0). If TRUE/1, 
#'  value of Ynode is given deterministically / constant.
#' @param f_gstar1 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  If a matrix/ data frame, its number of rows must be either nrow(data) or 1 (constant exposure assigned to all observations), and its number of 
#'  columns must be length(Anodes). If a vector, it must be of length nrow(data) or 1. If a function, it must return a data frame of counterfactual
#'  exposures sampled based on Anodes, Wnodes (and possibly Enodes and communityIndex) passed as a named argument "data". Thus, the function must 
#'  include "data" as one of its argument names. The interventions defined by f_gstar1 can be static, dynamic or stochastic. See Exmaples below.
#' @param f_gstar2 Either a function or a vector or a matrix/ data frame of counterfactual exposures, dependin on the number of exposure variables.
#'  It has the same components and requirements as f_gstar1
#' @param Qform Character vector of regression formula for Ynode. If not specified, the outcome variable is regressed on all covariates included in 
#'  Anodes, Wnodes and Enodes.
#' @param Qbounds Upper and lower bounds on Y and predicted values for initial Q. Defaults to the range of Y, widened by 10\% of the min and max values.
#' @param alpha Used to keep predicted values for initial Q bounded away from (0,1) for logistic fluctuation.
#' @param fluctuation Default to "logistic", it could also be "linear" (for targeting step).
#' @param f_g0 Optional function used to specify model knowledge about value of Anodes. It estimates \code{P(A | W, E)} under \code{g0} by 
#'  sampling a large vector/ data frame of Anode (of length or number of rows \code{nrow(data)*n_MCsims}) from \code{f_g0}
#' @param hform.g0 Character vector of regression formula for estimating the conditional density of P(A | W, E) under the observed treatment mechanism
#'  g0. If not specified, its form will be Anodes ~ Wnodes + Enodes. If there are more than one expsosure, it fits a joint probability.
#' @param hform.gstar Character vector of regression formula for estimating the conditional density P(A | W, E) under interventions f_gstar1 or f_gstar2. 
#'  If not specified, it follows the same rule used in hform.g0. 
#' @param lbound Value between (0,1) for truncation of predicted P(A | W, E). Default to 0.005
#' @param obs.wts Optional vector of observation (sampling) weights (of length nrow(data)). If NULL, assumed to be all 1. 
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

tmleCommunity <- function(data, Ynode, Anodes, Wnodes, Enodes = NULL, StratifyInd = NULL, YnodeDet = NULL, 
                          f_gstar1, f_gstar2 = NULL, Qform = NULL, Qbounds = NULL, alpha = 0.995, fluctuation = "logistic",                                                     
                          f.g0 = NULL, hform.g0 = NULL, hform.gstar = NULL, lbound = 0.005, obs.wts = NULL, 
                          h.g0_GenericModel = NULL, h.gstar_GenericModel = NULL, savetime.fit.hbars = TRUE, 
                          TMLE.targetStep = c("tmle.intercept", "tmle.covariate"),
                          n_MCsims = 1, 
                          CI_alpha = 0.05, 
                          rndseed = NULL, 
                          verbose = TRUE) {
  if (!is.null(rndseed))  set.seed(rndseed)  # make stochastic intervention trackable
  gvars$verbose <- verbose
  message("Running tmleCommunity with the following settings from tmleCom_Options(): "); str(gvars$opts)
  
  #----------------------------------------------------------------------------------
  # INITIALIZE PARAMETERS
  #----------------------------------------------------------------------------------
  if (is.null(savetime.fit.hbars)) savetime.fit.hbars <- getopt("savetime.fit.hbars")
  if (is.null(obs.wts)) obs.wts <- rep(1, nrow(data)) 
  if (!is.null(Qform) && !is.null(Ynode)) {
    Qform <- paste(Ynode, substring(Qform, first = as.numeric(gregexpr("~", Qform))))
    message("Since both Ynode and Qform are specified, the left-hand side of Qform will be ignored, with outcome being set to Ynode: " %+% Ynode)
    message("Thus the Qform becomes " %+% Qform)
  }
  if (!is.null(Qform) && is.null(Ynode)) {
    Ynode <- LhsVars(Qform)[1]
    message("Setting the Ynode to: " %+% Ynode)
  }
  TMLE.targetStep <- TMLE.targetStep[1]
  
  ## Check if any unexpected inputs
  if (!(TMLE.targetStep %in% c("tmle.intercept", "tmle.covariate"))) 
    stop("TMLE.targetStep argument must be either 'tmle.intercept' or 'tmle.covariate'")
  nodes <- list(Ynode = Ynode, Anodes = Anodes, Wnodes = Wnodes, Enodes = Enodes, Crossnodes = NULL)
  for (i in unlist(nodes)) {  CheckVarNameExists(data = data, varname = i) }
  if (!CheckInputs(data, nodes, Qform, hform.g0, hform.gstar, fluctuation, Qbounds, obs.wts)) stop()
  maptoYstar <- fluctuation=="logistic"  # if TRUE, cont Y values shifted & scaled to fall b/t (0,1)
  
  #----------------------------------------------------------------------------------
  # DEFINING (OPTIONAL) REGRESSION FORMS 
  #----------------------------------------------------------------------------------
  Q.sVars <- define_regform(as.formula(Qform), Anodes.lst = nodes$Ynode, Wnodes.lst = nodes[c("Anodes", "Wnodes", "Enodes")])
  h.g0.sVars <- define_regform(as.formula(hform.g0), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[3:5])
  if (!is.null(hform.gstar)) {
    h.gstar.sVars <- define_regform(as.formula(hform.gstar), Anodes.lst = nodes$Anodes, Wnodes.lst = nodes[c("Wnodes", "Enodes")])
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
  
  ## Create data based on Qform, hform.g0 and hform.gstar, in case of interaction or higher-order term.
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
  if (length(unique(obs.wts)) > 1 && any(unlist(OData$type.sVar[Anodes]) != "binary")) {
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
    TMLE.targetStep = TMLE.targetStep,
    obs.wts = obs.wts, 
    lbound = lbound,
    merged.form = merged.form, 
    model.Q.init = model.Q.init,
    Q.sVars = Q.sVars,
    f.g0 = f.g0,
    h.g0.sVars = h.g0.sVars,
    h.gstar.sVars = h.gstar.sVars,
    h.g0_GenericModel = h.g0_GenericModel,
    h.gstar_GenericModel = h.gstar_GenericModel,
    savetime.fit.hbars = savetime.fit.hbars,
    n_MCsims = n_MCsims
  ) 
  estinfo_list_g1 <- append(estinfo_list, list(f.gstar = f_gstar1))
  if (!is.null(f_gstar2)) { estinfo_list_g2 <- append(estinfo_list, list(f.gstar = f_gstar2))}
  
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
  EY_gstar1 <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, tmle_g_out = tmle_gstar1_out)
  EY_gstar2 <- NULL
  ATE <- NULL	
  otherInfo2 <- NULL
  if (!is.null(f_gstar2)) {
    EY_gstar2 <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, tmle_g_out = tmle_gstar2_out)
    ATE <- calcParameters(OData.ObsP0 = OData.ObsP0, inputYs = inputYs, alpha = CI_alpha, 
                          tmle_g_out = tmle_gstar1_out, tmle_g2_out = tmle_gstar2_out)
  }
  message("######################################################################################")
  message("Warning: inference for gcomp is not accurate! It is based on TMLE influence curves.")
  message("######################################################################################")
  
  tmleCommunity.res <- list(EY_gstar1 = EY_gstar1, EY_gstar2 = EY_gstar2, ATE = ATE)#, otherInfo1 = otherInfo1, otherInfo2 = otherInfo2)
  class(tmleCommunity.res) <- c(class(tmleCommunity.res), "tmleCommunity")
  return(tmleCommunity.res)
}
