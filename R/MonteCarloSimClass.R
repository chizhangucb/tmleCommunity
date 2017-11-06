#---------------------------------------------------------------------------------
# Monte-Carlo simulation sampling from stochastic g^* and evaluting G-Comp & TMLEs:
#---------------------------------------------------------------------------------


#------------------------------------ CalcMonteCarloEsts -----------------------------------
# Purpose: Evaluate G-Comp & TMLE using Monte-Carlo simulation. under stochastic intervention g^*
#-------------------------------------------------------------------------------------------
CalcMonteCarloEsts <- function(OData.ObsP0, OData.gstar, MC_fit_params, model.h.fit) {
  TMLE.targetStep <- MC_fit_params$TMLE.targetStep
  obs.wts <- MC_fit_params$obs.wts
  community.wts <- MC_fit_params$community.wts
  community.step <- MC_fit_params$community.step
  pooled.Q <- MC_fit_params$pooled.Q
  model.Q.init <- MC_fit_params$model.Q.init
  model.Q.star <- MC_fit_params$model.Q.star
  if (community.step != "NoCommunity") communityID <- OData.ObsP0$get.sVar(name.sVar = MC_fit_params$communityID)
  
  if (gvars$verbose) {
    message("================================================================")
    message("Running Monte Carlo evaluation of the substitution estimators...")
    message("================================================================")
  }
  evaluator <- MonteCarloSimClass$new(OData.ObsP0 = OData.ObsP0, OData.gstar = OData.gstar)
  nobs <- evaluator$nobs
  genMC.reps <- function(nrep)  {
    ## G-Comp estimator (i.e., GCOMP = MLE)
    MLE <- evaluator$get.gcomp(model.Q.init) # QY.init (G-Comp estimator) - est probY based on model for Q_Y
    ## TMLE estimator
    if (TMLE.targetStep == "tmle.covariate") {
      # TMLE A (adjusted by coefficient epsilon on h_bar ratio)
      TMLE <- evaluator$get.tmleCov(model.Q.star.cov = model.Q.star, model.h.fit = model.h.fit) 
    } else if (TMLE.targetStep == "tmle.intercept") {  
      # TMLE B (adjusted by intercept epsilon where h_bar were used as weights)
      TMLE <- evaluator$get.tmleInt(model.Q.star.int = model.Q.star) 
    } else {
      TMLE <- NA
    }
    ## fi_W - hold W fixed to observed values (a component of TMLE Var)
    fiWs_list <- evaluator$get.fiW()  
    # Put all estimators together and add names (defined in out_nms outside of this function):
    if (community.step == "individual_level") {
      if (!is.null(communityID)) {
        TMLE <- aggregate(x = TMLE, by=list(newid = communityID), mean)
        MLE <- aggregate(x = MLE, by=list(newid = communityID), mean)[, 2]
        obs.wts <- community.wts[match(TMLE[, "newid"], community.wts[, "id"]), "weights"]
        TMLE <- TMLE[, 2]
      } else {
        warningMesg <- c("Since individual-level TMLE requires communityID to aggregate data to the cluster-level in the end. ",
                         "Lack of 'communityID' forces the algorithm to treat the data as non-hierarchical.")
        warning(warningMesg[1] %+% warningMesg[2])
      }
    }
    mean_psis_all <- c(TMLE = weighted.mean(TMLE, w = obs.wts), 
                       MLE = weighted.mean(MLE, w = obs.wts), 
                       fiWs_list$fiW_Qinit)
    names(mean_psis_all) <- c("TMLE", "MLE", paste("fWi_init_", c(1:nobs), sep = ""))
    return(mean_psis_all)
  }
  psi_est_mean <- genMC.reps(1)
  return(psi_est_mean)
}


## ---------------------------------------------------------------------
#' R6 class for evaluating different plug-in estimators via Monte-Carlo resampling where new exposures are generated 
#' under the user-specified arbitrary intervention function g*.
#' 
#'  \code{MonteCarloSimClass} only resamples \code{A} under the intervention function \code{g_star}, never \code{W} or \code{E}.
#'  For each MC simulation, it firstly treats \code{model.Q.init} as the fitted models for \code{E[Y|A,W,E]}, then estimate 
#'  \code{psi_n} using Monte-Carlo integration. i.e., average of \code{n} predicted \code{E(Y|A=a, W=w,E=e)} where \code{a} is 
#'  a vector of \code{n} new exposures randomly drawn under \code{g_star}. Take as many iterations as needed until convergence 
#'  of \eqn{\psi^{I}_n} (or \eqn{\psi^{II}_n}) occurs. 
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{OData.ObsP0}} - A \code{DatKeepClass} class object, where exposures are generated under observed exposure mechanism g0.
#' \item{\code{OData.gstar}} - A \code{DatKeepClass} class object, where new exposures are generated under user-specified 
#'  intervention \eqn{g^{*}}.
#' \item{\code{model.Q.init}} - A fitted model for \code{E[Y|A,W,E]}.
#' \item{\code{model.Q.star.cov}} - A targeting model for covariate-based unweighted TMLE.
#' \item{\code{model.Q.star.int}} - A targeting model for intercept-based TMLE.
#' \item{\code{QY.init}} - Estimates of G-COMP.
#' \item{\code{QY.star.cov}} - Estimates of covariate-based unweighted TMLE.
#' \item{\code{QY.star.int}} - Estimates of intercept-based TMLE.
#' \item{\code{nobs}} - Number of observations in the observed data frame.
#' \item{\code{p}} - Number of Monte-Carlo simulations performed.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(OData.ObsP0, OData.gstar, ...)}}{Instantiate an new instance of \code{MonteCarloSimClass}.}
#'   \item{\code{get.gcomp(m.Q.init)}}{Predict \code{QY.init} = \eqn{E[Y_{g^*}]} based on the initial model fit \code{model.Q.init}.}
#'   \item{\code{get.tmleCov(model.Q.star.cov, model.h.fit)}}{Update \code{QY.init} based on the targeting model \code{model.Q.star.cov}
#'    and the model for clever covriate h \code{model.h.fit}.}
#'   \item{\code{get.tmleCov(model.Q.star.cov, model.h.fit)}}{Update \code{QY.init} based on the targeting model \code{model.Q.star.cov}
#'    and the model for clever covriate h \code{model.h.fit}.}
#'   \item{\code{get.tmleInt(model.Q.star.int)}}{Update \code{QY.init} based on the targeting model \code{model.Q.star.int}.}
#'   \item{\code{get.fiW()}}{Get an estimate of fiW (hold ALL W's fixed).}
#' }
#' @import data.table
#' @export
MonteCarloSimClass <- R6Class(classname = "MonteCarloSimClass",
  portable = TRUE,
  class = TRUE,
  public = list(
    OData.ObsP0 = NULL, 
    OData.gstar = NULL,
    model.Q.init = NULL,
    model.Q.star.cov = NULL,
    model.Q.star.int = NULL,
    QY.init = NULL,
    QY.star.cov = NULL,
    QY.star.int = NULL,
    nobs = NA_integer_,        # no. of samples in the OBSERVED (original) data
    p = NA_integer_,           # no. of times n-size A were resampled under gstar
    
    initialize = function(OData.ObsP0, OData.gstar, ...) {
      self$OData.ObsP0 <- OData.ObsP0
      self$nobs <- OData.ObsP0$nobs
      self$OData.gstar <- OData.gstar
      self$p <- OData.gstar$p
      invisible(self)
    },
    
    # MLE - Predict E[Yg_star] (QY.init) for each i, based on the initial model fit for E[Y|C^Y] (model.Q.init)
    # output is a vector of length n*p
    get.gcomp = function(model.Q.init) {
      self$model.Q.init <- model.Q.init
      OData.ObsP0 <- self$OData.ObsP0
      OData.gstar <- self$OData.gstar
      # print("getting predictions for gcomp...")
      # predictions P(Y=1) under A^* for non-det Y in OData.ObsP0
      QY.init <- model.Q.init$predict(newdata = OData.gstar)$getprobA1  
      QY.init[OData.ObsP0$det.Y] <- OData.ObsP0$noNA.Ynodevals[OData.ObsP0$det.Y]
      self$QY.init <- QY.init
      invisible(QY.init)
    },
    
    # TMLE Covariate - Update QY.init based on the est. coefficient for the clever covariate h_g0/h_gstar in univar. model fluct
    # output is a vector of length n*p
    get.tmleCov = function(model.Q.star.cov, model.h.fit) {
      if (is.null(self$QY.init)) stop("call MonteCarloSimClass$get.gcomp first") # QY.init <- self$get.gcomp(self$model.Q.init)
      iptw <- predict.hbars(newdatnet = self$OData.gstar, model.h.fit = model.h.fit)
      if (!is.na(coef(model.Q.star.cov))) {
        self$model.Q.star.cov <- model.Q.star.cov
        off <- qlogis(self$QY.init)
        self$QY.star.cov <- plogis(off + coef(model.Q.star.cov)*iptw)
      } else {
        self$QY.star.cov <- self$QY.init
      }
      invisible(self$QY.star.cov)
    },
    
    # TMLE Intercept - Update QY.init based on the est. intercept of the model fluct (m.Q.star_h_B)
    # output is a vector of length n*p
    get.tmleInt = function(model.Q.star.int) {
      if (is.null(self$QY.init)) stop("call MonteCarloSimClass$get.gcomp first") # QY.init <- self$get.gcomp(self$model.Q.init)
      if (!is.na(coef(model.Q.star.int))) {
        self$model.Q.star.int <- model.Q.star.int
        off <- qlogis(self$QY.init)
        self$QY.star.int <- plogis(off + coef(model.Q.star.int))
      } else {
        self$QY.star.int <- self$QY.init
      }
      invisible(self$QY.star.int)
    },
    
    # Get an estimate of fiW (hold ALL W's fixed at once) - a component of TMLE Var
    # output is a vector of length n*p, where each of n obs is then averaged p times.
    get.fiW = function() {
      if (is.null(self$QY.init)) stop("call MonteCarloSimClass$get.gcomp first") # QY.init <- self$get.gcomp(self$m.Q.init)
      # ******* fi_W based on Q, N.init model ******
      ID <- rep.int(c(1 : self$nobs), self$p)
      # taking average over p samples for each of n obs
      # data.table is used on fiW_Qinit for performing unit-level mean and var evaluation (n-length result)
      fiW_Qinit <- data.table::data.table(ID = ID, fiW = self$QY.init)  # Use P(Y=1) under A^* for non-det Y in OData.ObsP0
      fiW_Qinit.mean <- fiW_Qinit[, lapply(.SD, mean, na.rm=TRUE), by="ID", .SDcols=c("fiW") ][["fiW"]]
      fiW_Qinit.var <- fiW_Qinit[, lapply(.SD, var, na.rm=TRUE), by="ID", .SDcols=c("fiW") ][["fiW"]]
      # If p=1, no mean or var performed and fiW_Qinit.mean = self$QY.init (Results from get.gcomp) & fiW_Qinit.var = NA
      return(list(fiW_Qinit = fiW_Qinit.mean, fiW_Qinit.var = fiW_Qinit.var))
    }
  ),
    
  active = list(
    # plchld = function() {}
  ),
    
  private = list(
    # placeholder = list()
  )
)
