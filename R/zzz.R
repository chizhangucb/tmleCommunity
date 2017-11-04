#-----------------------------------------------------------------------------
# Global State Vars (can be controlled globally with options(stremr.optname = ))
#-----------------------------------------------------------------------------
gvars <- new.env(parent = emptyenv())
gvars$verbose <- FALSE      # verbose mode (print all messages)
gvars$opts <- list()        # named list of package options that is controllable by the user (set_all_stremr_options())
gvars$misval <- NA_integer_ # the default missing value for observations (# gvars$misval <- -.Machine$integer.max)
gvars$misXreplace <- 0L     # the default replacement value for misval that appear in the design matrix
gvars$tolerr <- 10^-12      # tolerance error: assume for abs(a-b) < gvars$tolerr => a = b
gvars$sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")

getopt <- function(optname) {
  opt <- gvars$opts
  if (!(optname %in% (names(opt)))) stop(optname %+% ": this options does not exist")
  return(opt[[optname]])
}

`%+%` <- function(a, b) paste0(a, b)
  
gvars$opts.allowedVals <- list(Qestimator = c("speedglm__glm", "glm__glm", "h2o__ensemble", "SuperLearner"),
                               gestimator = c("speedglm__glm", "glm__glm", "h2o__ensemble", "SuperLearner"),
                               # fitclass = c("speedglmS3", "glmS3", "h2oS3", "SLS3"),
                               bin.method = c("equal.mass", "equal.len", "dhist"),
                               nbins = "_positive_integer_",
                               maxncats = "_positive_integer_",  # Max number of unique categories a categorical variable A[j] can have.  
                                                                 # If A[j] has more, it is automatically considered as continuous.
                               maxNperBin = "_positive_integer_",
                               parfit = c(TRUE, FALSE),
                               poolContinVar = c(TRUE, FALSE),
                               savetime.fit.hbars = c(TRUE, FALSE),
                               h2ometalearner = "_character_",
                               h2olearner = "_character_",
                               CVfolds = "_positive_integer_",
                               SL.library = "_character_"
                               # panel.effect = c("individual", "time", "twoways"), 
                               # panel.model = c("within", "random", "between","pooling", "ht", "fd"),
                               # community.index = "_character_"
  )

#' Print Current Option Settings for \code{tmleCommunity}
#'
#' Print Current Option Settings for \code{tmleCommunity}
#' @return Invisibly returns a list of \code{tmleCommunity} options.
#' @seealso \code{\link{tmleCom_Options}}
#' @export
print_tmleCom_opts <- function() {
  print(gvars$opts)
  invisible(gvars$opts)
}

#' Setting all possible \code{tmleCommunity} options
#'
#' Additional options that control \code{tmleCommunity} package
#' @param Qestimator A string specifying default estimator for outcome mechanism model fitting. 
#'  The default estimator is \code{"speedglm__glm"}, which estimates regressions with \code{\link[speedglm]{speedglm.wfit}}; 
#'  Estimator \code{"glm__glm"} uses \code{\link[stats]{glm.fit}};
#'  Estimator \code{"h2o__ensemble"} implements the super learner ensemble (stacking) algorithm using the H2O R interface; 
#'  Estimator \code{"SuperLearner"} implements the super learner prediction methods.
#'  Note that if \code{"h2o__ensemble"} fails, it falls back on {"SuperLearner"}. If \code{"SuperLearner"} fails, 
#'  it falls back on {"speedglm__glm"}. If \code{"speedglm__glm"} fails, it falls back on {"glm__glm"}.
#' @param gestimator A string specifying default estimator for exposure mechanism fitting. It has the same options as \code{Qestimator}.
#' @param bin.method Specify the method for choosing bins when discretizing the conditional continuous exposure variable \code{A}.
#'  The default method is \code{"equal.mass"}, which provides a data-adaptive selection of the bins based on equal mass/ area, i.e., 
#'  each bin will contain approximately the same number of observations as otheres. Method \code{"equal.len"} partitions the range of 
#'  \code{A} into equal length \code{nbins} intervals. Method \code{"dhist"} uses a combination of the above two approaches. Please
#'  see Denby and Mallows "Variations on the Histogram" (2009) for more details. Note that argument \code{maxNperBin} controls
#'  the maximum number of observations in each bin.
#' @param nbins When \code{bin.method = "equal.len"}, set to the user-supplied number of bins when discretizing a continous variable/
#'  If not specified, then default to 5; If setting to as \code{NA}, then set to the nearest integer of \code{nobs/ maxNperBin}, where
#'  \code{nobs} is the total number of observations in the input data. When method is \code{"equal.mass"}, \code{nbins} will be set as 
#'  the maximum of the default \code{nbins} and the nearest integer of \code{nobs/ maxNperBin}.
#' @param maxncats Integer that specifies the maximum number of unique categories a categorical variable \code{A[j]} can have. If 
#'  \code{A[j]} has more unique categories, it is automatically considered a continuous variable. Default to 10.
#' @param maxNperBin Integer that specifies the maximum number of observations in each bin when discretizing a continuous variable 
#'  \code{A[j]} (applies directly when \code{bin.method =} \code{"equal.mass"} and indirectly when \code{bin.method = "equal.len"}, but 
#'  \code{nbins = NA}).
#' @param parfit Logical. If \code{TRUE}, perform parallel regression fits and predictions for discretized continuous variables by 
#'  functions \code{\link{foreach}} and \code{\link{dopar}} in \code{foreach} package. Default to \code{FALSE}. Note that it requires 
#'  registering a parallel backend prior to running \code{tmleCommunity} function, e.g., using \code{doParallel} R package and running 
#'  \code{registerDoParallel(cores = ncores)} for \code{ncores} parallel jobs.
#' @param poolContinVar Logical. If \code{TRUE}, when fitting a model for binirized continuous variable, pool bin indicators across
#'  all bins and fit one pooled regression. Default to \code{FALSE}.
#' @param savetime.fit.hbars Logical. If \code{TRUE}, skip estimation and prediction of exposure mechanism P(A|W,E) under g0 & gstar
#'  when \code{f.gstar = NULL} and \code{TMLE.targetStep = "tmle.intercept"}, and then directly set \code{h_gstar_h_gN = 1} for each 
#'  observation. Default to \code{TRUE}.
#' @param h2ometalearner A string to pass to \code{\link{h2o.ensemble}}, specifying the prediction algorithm used to learn the optimal 
#'  combination of the base learners. Supports both h2o and SuperLearner wrapper functions. Default to "h2o.glm.wrapper".  
#' @param h2olearner A string or character vector to pass to \code{\link{h2o.ensemble}}, naming the prediction algorithm(s) used to train
#'  the base models for the ensemble. The functions must have the same format as the h2o wrapper functions. Default to "h2o.glm.wrapper".
#' @param CVfolds Set the number of splits for the V-fold cross-validation step to pass to \code{\link{SuperLearner}} and 
#'  \code{\link{h2o.ensemble}}. Default to 5.
#' @param SL.library A string or character vector of prediction algorithms to pass to \code{\link{SuperLearner}}. Default to 
#'  c("SL.glm", "SL.step", "SL.glm.interaction"). For more available algorithms see \code{SuperLearner::listWrappers()} .
#' @return Invisibly returns a list with old option settings.
#' 
#' @seealso \code{\link{print_tmlenet_opts}}
#' 
#' @examples
#' tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", nbins = 20)
#'
#' tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = 10000,
#'                     h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), CVfolds = 10)
#'
#' tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
#'                     SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5)
#' 
#' tmleCom_Options(Qestimator = "SuperLearner", gestimator = "h2o__ensemble", maxNperBin = nrow(data),
#'                     SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
#'                 h2ometalearner = "h2o.deeplearning.wrapper", 
#'                 h2olearner = c("h2o.gbm.wrapper", "h2o.randomForest.wrapper"))
#' 
#' @export
tmleCom_Options <- function(Qestimator = c("speedglm__glm", "glm__glm", "h2o__ensemble", "SuperLearner"),  
                            gestimator = c("speedglm__glm", "glm__glm", "h2o__ensemble", "SuperLearner"),  
                            bin.method = c("equal.mass", "equal.len", "dhist"),
                            nbins = 5, 
                            maxncats = 10,  
                            maxNperBin = 500,
                            parfit = FALSE,
                            poolContinVar = FALSE,
                            savetime.fit.hbars = TRUE,
                            h2ometalearner = "h2o.glm.wrapper",
                            h2olearner = "h2o.glm.wrapper",
                            CVfolds = 5,
                            SL.library = c("SL.glm", "SL.step", "SL.glm.interaction")
                           ) {
  old.opts <- gvars$opts
  Qestimator <- Qestimator[1L]
  gestimator <- gestimator[1L]
  # fitclass <- fitclass[1L]
  bin.method <- bin.method[1L]
  # panel.effect <- panel.effect[1L]
  # panel.model <- panel.model[1L]

  if (!(Qestimator %in% gvars$opts.allowedVals[["Qestimator"]])) 
    stop("Qestimator must be one of: " %+% paste0(gvars$opts.allowedVals[["Qestimator"]], collapse=", "))
  if (!(gestimator %in% gvars$opts.allowedVals[["gestimator"]])) 
    stop("gestimator must be one of: " %+% paste0(gvars$opts.allowedVals[["gestimator"]], collapse=", "))
  if (!(bin.method %in% gvars$opts.allowedVals[["bin.method"]])) 
    stop("bin.method must be one of: " %+% paste0(gvars$opts.allowedVals[["bin.method"]], collapse=", "))
  
  if (any(c(Qestimator, gestimator) %in% "h2o__ensemble")) {
    if (!requireNamespace("h2o") || !requireNamespace("h2oEnsemble")) 
      stop("h2o and h2oEnsemble package are required if either Qestimator or gestimator is 'h2o__ensemble'")
  }
  if (any(c(Qestimator, gestimator) %in% "SuperLearner")) {
    if (!requireNamespace("SuperLearner"))  stop("SuperLearner package is required if either Qestimator or gestimator is 'SuperLearner'.")
  }

  # if (Qestimator == "speedglm__glm") { Qfitclass <- "speedglmS3" }
  # if (Qestimator == "glm__glm") { Qfitclass <- "glmS3" }
  # if (Qestimator == "h2o__ensemble") { Qfitclass <- "h2oS3" }  
  # if (Qestimator == "SuperLearner") { Qfitclass <- "SLS3" }  
    
  opts <- list(
    Qestimator = Qestimator,
    gestimator = gestimator,
    bin.method = bin.method,
    nbins = nbins,
    maxncats = maxncats,
    maxNperBin = maxNperBin,
    parfit = parfit,
    poolContinVar = poolContinVar,
    savetime.fit.hbars = savetime.fit.hbars,
    h2ometalearner = h2ometalearner,
    h2olearner = h2olearner,
    CVfolds = CVfolds,
    SL.library = SL.library
  )
  gvars$opts <- opts
  invisible(old.opts)
}

# returns a function (alternatively a call) that tests for missing values in (sA, sW)
testmisfun <- function() {
  if (is.na(gvars$misval)) {
    return(is.na)
  } else if (is.null(gvars$misval)){
    return(is.null)
  } else if (is.integer(gvars$misval)) {
    return(function(x) {x==gvars$misval})
  } else {
    return(function(x) {x%in%gvars$misval})
  }
}
    
get.misval <- function() {
  gvars$misfun <- testmisfun()
  gvars$misval
}
    
set.misval <- function(gvars, newmisval) {
  oldmisval <- gvars$misval
  gvars$misval <- newmisval
  gvars$misfun <- testmisfun()    # EVERYTIME gvars$misval HAS CHANGED THIS NEEDS TO BE RESET/RERUN.
  invisible(oldmisval)
}
gvars$misfun <- testmisfun()

# Allows tmleCommunity functions to use e.g., getOption("tmleCommunity.verbose") to get verbose printing status
#.onLoad <- function(libname, pkgname) {
#  op <- options()
#  op.tmleCommunity <- list(
#    tmleCommunity.verbose = gvars$verbose
#  )
#  # reset all options to their defaults on load:
#  tmleCommunity_options()
#  toset <- !(names(op.tmleCommunity) %in% names(op))
#  if(any(toset)) options(op.tmleCommunity[toset])
#  invisible()
#}
