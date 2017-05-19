# ---------------------------------------------------------------------------------
# Generic S3 constructors for the summary model classes:
# ---------------------------------------------------------------------------------
newsummarymodel <- function(reg, DatKeepClass.g0, ...) { UseMethod("newsummarymodel") }
# Summary model constructor for generic regression with multivariate outcome, but one set of predictors
newsummarymodel.generic <- function(reg, DatKeepClass.g0, ...) GenericModel$new(reg = reg, DatKeepClass.g0 = DatKeepClass.g0, ...)
# Summary model constructor for continuous outcome A[j]:
newsummarymodel.contin <- function(reg, DatKeepClass.g0, ...) ContinModel$new(reg = reg, DatKeepClass.g0 = DatKeepClass.g0, ...)
# Summary model constructor for categorical outcome A[j]:
newsummarymodel.categor <- function(reg, DatKeepClass.g0, ...) CategorModel$new(reg = reg, DatKeepClass.g0 = DatKeepClass.g0, ...)
# Summary model constructor for binary outcome A[j]:
newsummarymodel.binary <- function(reg, ...) BinaryOutModel$new(reg = reg, ...)


## ---------------------------------------------------------------------
#' R6 class for defining regression models that evaluate multivariate joint conditional density P(A|W,E) (or P(A|W) if community-specific)
#'
#'  \code{RegressionClass} provides multiple options used when estimating a joint conditional density \code{P(A|W,E)}. Note that \code{A} 
#'  can be multivariate, if so, hazard specification will factorize \code{P(A|W,E)} = \code{P(A[1],...,A[j]|W,E)} as a sequence
#'  \code{P(A[1]|W,E)} * \code{P(A[2]|W, E, A[1])} * ... *               
#'  \code{P(A[j]|W, E, A[1],...,A[j-1])}, where each of the compoenents \code{A[i]} can be 
#'  either binary, categorical or continuous, and each of the conditional densities \code{P(A[i]|W, E, A[1],...,A[i-1])} will be controlled by 
#'  a new instance of \code{\link{GenericModel}} class. If \code{A[i]} binary, \code{P(A[i]|W, E, A[1],...,A[i-1])} will be esimtated by a user-
#'  specific library of candidate algorithms, including parametric estimators such as logistic model with only main terms, and data-adaptive 
#'  estimator such as super-learner algorithms. If \code{A[i]} continuous (or categorical), \code{P(A[i]|W, E, A[1],...,A[i-1])} will then be 
#'  controlled by a new instance of \code{\link{ContinModel}} class (or \code{\link{CategorModel}} class). Note that each \code{GenericModel}, 
#'  \code{ContinModel} and \code{CategorModel} class will accompany with an adjunctive clone of a parent \code{RegressionClass} class. The  
#'  automatically recursive process of defining new instances of \code{GenericModel} and cloning \code{RegressionClass} classes won't stop until
#'  the entire sequence of binary regressions that represents \code{P(A|W,E)} is constructed.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{outvar.class}} - Character vector of outcome variable classes (of length(outvar)): one of \code{bin}, \code{cont}, \code{cat}.
#' \item{\code{outvar}} - Character vector of regression outcome variable names.
#' \item{\code{predvars}} - Character vector of regression-specific predictor names or a pool of all available predictor names.
#' \item{{reg_hazard}} - Logical, hazard fitting method. If TRUE, factorize P(outvar | predvars) into \prod_{j}{P(outvar[j] | predvars)} for each j.
#' \item{\code{subset_vars}} - Named list for subset variables/ expression (later converted to logical vector).
#' \item{\code{ReplMisVal0}} - Logical. If TRUE, replace all gvars$misval among predictors with user-supplied gvars$misXreplace (Default to 0).
#' \item{\code{nbins}} - Number of bins used for estimation of a continuous outvar, defined in ContinModel$new().
#' \item{\code{estimator}} - Character, one of "speedglm__glm" (default), "glm__glm", "h2o__ensemble", "SuperLearner". The estimator for which to fit 
#'    regression model. For "h2o__ensemble" and "SuperLearner", users can specify the data-adaptive algorithms through \code{tmleCom_Options}.
#' \item{\code{parfit}} - Logical. If TRUE, use parallel computing on binary regressions. See \code{foreach::foreach}.
#' \item{\code{pool_cont}} - Logical. If TRUE, pool over bins of a continuous outvar and fit one regression, along with bin_ID as an extra variable. 
#' \item{\code{outvars_to_pool}} - Character vector of bin names of a continuous outvars, should be identical to \code{bin_nms}.
#' \item{\code{intrvls}} - Numeric vector defining the number and positions of the bins or a named list of numeric vectors if 2 or more outvars.
#'    If not specified and outvar continuous, intervals will be determined in \code{ContinModel} through \code{DatKeepClass$detect.sVar.intrvls}.
#' \item{\code{intrvls.width}} - Named numeric vector of bin widths for each bin in \code{self$intrvls}. If not specified, default to 1 if
#'    outvar binary, default to \code{diff(self$intrvls)} if outvar continuous, 
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(outvar.class = gvars$sVartypes$bin,
#'                   outvar, predvars, subset_vars, intrvls,
#'                   ReplMisVal0 = TRUE,
#'                   estimator = getopt("Qestimator"),
#'                   parfit = getopt("parfit"),
#'                   pool_cont = getopt("poolContinVar")}}{Instantiate an new instance of \code{RegressionClass}}
#'   \item{\code{ChangeManyToOneRegresssion(k_i, reg)}}{Clone the parent \code{RegressionClass} (\code{reg}) that include \code{length(self$outvar)} 
#'   regressions, and reset self to a single univariate \code{k_i} regression for outcome \code{self$outvar[[k_i]]}.}
#'   \item{\code{resetS3class()}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{S3class(newclass)}}{...}
#'   \item{\code{get.reg}}{...}
#' }
#' @export
RegressionClass <- R6Class("RegressionClass",
  class = TRUE,
  portable = TRUE,
  public = list(
    outvar.class = character(),    # Vector for classes of the outcome vars: bin / cont / cat
    outvar = character(),          # Vector of regression outcome variable names
    predvars = character(),        # Either a pool of all predictors (W) or regression-specific predictor names
    reg_hazard = FALSE,            # If TRUE, the joint P(outvar|predvars) is factorized as \prod_{j}{P(outvar[j] | predvars)} for each j outvar 
    subset_vars = NULL,            # Named LIST for subset vars (later evaluated to logical vector in the envir of the data)
    ReplMisVal0 = TRUE,            # If TRUE all gvars$misval among predicators are replaced with with gvars$misXreplace (0) 
    nbins = NULL,                  # Actual nbins used, for cont. outvar, defined in ContinModel$new()
    estimator = character(),       # Specify default estimator for model fitting
    parfit = logical(),            # TRUE for fitting binary regressions in parallel
    pool_cont = logical(),         # Pool binned cont outvar obs into long format (adding bin_ID as a covaraite)
    outvars_to_pool = character(), # Names of the binned continuous sVars, should match bin_nms
    intrvls = numeric(),           # Vector of numeric cutoffs defining the bins or a named list of numeric intervals (for length(self$outvar) > 1)
    intrvls.width = 1L,            # Named vector of bin-widths (bw_j : j=1,...,M) for each each bin in self$intrvls
	                           # When A is not continuous, intrvls.width IS SET TO 1.
                                   # When A is continuous, intrvls.width is SET TO self$intrvls.width 
    initialize = function(outvar.class = gvars$sVartypes$bin,
                          outvar, predvars, subset_vars, intrvls,
                          ReplMisVal0 = TRUE, # add ReplMisVal0 = TRUE for case sA = (netA, sA[j]) with sA[j] continuous, causing an error otherwise:
                          estimator = getopt("Qestimator"),  # By defautl, g and Q use the same estimator
                          parfit = getopt("parfit"),
			  nbins = getopt("nbins"),
                          pool_cont = getopt("poolContinVar")
                          ) {
      assert_that(length(outvar.class) == length(outvar))

      self$outvar.class <- outvar.class
      self$outvar <- outvar
      self$predvars <- predvars
      self$ReplMisVal0 <- ReplMisVal0
      self$estimator <- estimator
      self$parfit <- parfit
      self$nbins <- nbins
      self$pool_cont <- pool_cont

      n_regs <- length(outvar)

      if (!missing(intrvls)) {
        assert_that(is.list(intrvls))
        assert_that(length(outvar) == length(intrvls))
        assert_that(all(names(intrvls) %in% outvar))
        self$intrvls <- intrvls
      } else {
        self$intrvls <- NULL
      }

      if (!missing(subset_vars)) {
        self$subset_vars <- subset_vars
        if (length(subset_vars) < n_regs) {
          self$subset_vars <- rep_len(subset_vars, n_regs)
        } else if (length(subset_vars) > n_regs) {
          if (!is.logical(subset_vars)) stop("not implemented")
	  self$subset_vars <- subset_vars[1:n_regs]  # chose the first n_regs elements in subset_vars
          # ... TO FINISH ... increase n_regs to all combinations of (n_regs x subset_vars)
        }
      } else {
        self$subset_vars <- rep_len(list(TRUE), n_regs)
      }
    },

    # take the clone of a parent RegressionClass (reg) for length(self$outvar) regressions
    # and set self to a single univariate k_i regression for outcome self$outvar[[k_i]]
    ChangeManyToOneRegresssion = function(k_i, reg) {
      assert_that(!missing(k_i))
      if (missing(reg)) stop("reg must be also specified when k_i is specified")
      assert_that(is.count(k_i))
      assert_that(k_i <= length(reg$outvar))

      n_regs <- length(reg$outvar)
      self$outvar.class <- reg$outvar.class[[k_i]] # Class of the outcome var: binary, categorical, continuous:
      self$outvar <- reg$outvar[[k_i]] # An outcome variable that is being modeled:

      if (self$reg_hazard) {
        # when modeling bin hazard indicators, no need to condition on previous outcomes as they will all be degenerate
        self$predvars <- reg$predvars # Regression covars (predictors):
      } else {
        self$predvars <- c(reg$outvar[-c(k_i:n_regs)], reg$predvars) # Regression covars (predictors):
      }

      # the subset is a list when RegressionClass specifies several regression models at once,
      # obtain the appropriate subset for this regression k_i and set it to self
      if (is.list(reg$subset_vars)) {
        self$subset_vars <- reg$subset_vars[[k_i]]
      }
      if (is.list(reg$intrvls)) {
        outvar_idx <- which(names(reg$intrvls) %in% self$outvar)
        self$intrvls <- reg$intrvls[[outvar_idx]]
      }
      self$S3class <- self$outvar.class # Set class on self for S3 dispatch...
      return(invisible(self))
    },

    resetS3class = function() class(self) <- c("RegressionClass", "R6")
  ),

  active = list(
    # For S3 dispatch on newsummarymodel():
    S3class = function(newclass) {
      if (!missing(newclass)) {
        if (length(newclass) > 1) stop("cannot set S3 class on RegressionClass with more than one outvar variable")
        if (length(class(self)) > 2) stop("S3 dispatch class on RegressionClass has already been set")
        class(self) <- c(class(self), newclass)
      } else {
        return(class(self))
      }
    },

    get.reg = function() {
      list(outvar.class = self$outvar.class,
          outvar = self$outvar,
          predvars = self$predvars,
          subset_vars = self$subset_vars)
    }
  )
)
	
	
## ---------------------------------------------------------------------
#' Generic R6 class for modeling (fitting and predicting) P(A=a | W=w, E=e) where A can be multivariate (A[1], ..., A[j])
#'
#'  \code{GenericModel} defines and models the conditional density \code{P(A=a|W=w,E=e)}, where \code{a} are generated under \code{g_star}
#'  or \code{g_0}. By calling \code{self$new(reg)}, it utilizes estimation options defined in \code{\link{RegressionClass}} class, and 
#'  automatically factorizes \code{P(A|W,E)} into an entire tree of binary regression models, where a new instance of \code{\link{BinaryOutModel}}
#'  class will be initialized for each binary regression. By calling \code{self$fit(data)} and \code{self$predict(newdata)}, where \code{data} and    
#'  \code{newdata} are \code{\link{DatKeepClass}} class objects, it fits \code{P(A|W,E)} and predicts \code{P(A=1|W=w,E=e)}, where values of 
#'  {\code{w}, \code{e}} are from \code{newdata}. Moreover, it predicts likelihood function P(A=a| W=w,E=e) through \code{self$predictAeqa(newdata)},
#'  where {\code{a}, \code{w}, \code{e}} are from \code{newdata} (also a \code{\link{DatKeepClass}} class).
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{n_regs}} - Total number of regression models.
#' \item{\code{parfit_allowed}} - Logical. If TRUE, allow parallel computing to fit multivariate outvar.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, ...)}}{Use \code{reg} (a \code{\link{RegressionClass}} class object) to instantiate an new object of \code{GenericModel}}
#'   \item{\code{length}}{...}
#'   \item{\code{getPsAsW.models}}{Get all model objects (one model object per outcome var A[j])}
#'   \item{\code{getcumprodAeqa}}{Get joint prob as a vector of the cumulative prod over j for P(A[j]=a[j]|W,E)}
#'   \item{\code{fit(data, savespace = TRUE)}}{savespace is Logical. If TRUE, wipe out all internal data when doing many stacked regressions}
#'   \item{\code{copy.fit(Generic.Model)}}{...}
#'   \item{\code{predict(newdata, savespace = TRUE)}}{...}
#'   \item{\code{predictAeqa(newdata, savespace = TRUE, wipeProb = TRUE, ...)}}{calculate the likelihood for observed A, W, E. Note that wipeProb is 
#'     logical argument for \code{self$wipe.alldat}. If FALSE, vectors of probA1 & probAeqa will be kept.}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{wipe.alldat(wipeProb = TRUE)}}{...}
#' }
#' @export
GenericModel <- R6Class(classname = "GenericModel",
  portable = TRUE,
  class = TRUE,
  public = list(
    n_regs = integer(),        # Total no. of reg. models (Nonparametric / logistic regressions)
    parfit_allowed = FALSE,    # Allow parallel fit of multivar outvar when 1) reg$parfit = TRUE & 2) all.outvar.bin = TRUE
    is.fitted = FALSE,         
    initialize = function(reg, ...) {
      self$n_regs <- length(reg$outvar) # Number of sep. logistic regressions to run
      all.outvar.bin <-  all(reg$outvar.class %in% gvars$sVartypes$bin)

      if (reg$parfit & all.outvar.bin & (self$n_regs > 1)) self$parfit_allowed <- TRUE

      if (gvars$verbose) {
        print("#----------------------------------------------------------------------------------")
        print("New instance of GenericModel:")
        print("#----------------------------------------------------------------------------------")
        print("Outcomes: " %+% paste(reg$outvar, collapse = ", "))
        print("Predictors: " %+% paste(reg$predvars, collapse = ", "))
        print("No. of regressions: " %+% self$n_regs)
        print("All outcomes binary? " %+% all.outvar.bin)
        if (self$parfit_allowed) print("Performing parallel fits: " %+% self$parfit_allowed)
        print("#----------------------------------------------------------------------------------")
      }

      # factorize the joint into univariate regressions, by dimensionality of the outcome variable (A_nms):
      for (k_i in 1:self$n_regs) {
	reg_i <- reg$clone()
	reg_i$ChangeManyToOneRegresssion(k_i, reg)
        # Calling the constructor for the summary model P(A[j]|\bar{A}[j-1], W}), dispatching on reg_i class
	regS3class <- reg_i$S3class
        if (is.null(regS3class)) {
          reg_i$S3class <- "generic"; class(reg_i$S3class) <- "generic"
        }
        PsAsW.model <- newsummarymodel(reg = reg_i, ...)
        private$PsAsW.models <- append(private$PsAsW.models, list(PsAsW.model))
        names(private$PsAsW.models)[k_i] <- "P(A|W)."%+%k_i
      }
      invisible(self)
    },
    length = function(){ base::length(private$PsAsW.models) },
    getPsAsW.models = function() { private$PsAsW.models },  # get all summary model objects (one model object per outcome var A[j])
    getcumprodAeqa = function() { private$cumprodAeqa },  # get joint prob as a vector of the cumulative prod over j for P(A[j]=a[j]|W)

    # -------------------------------------------------------------------------------------
    # Methods for fitting regression models - fit, copy.fit
    # -------------------------------------------------------------------------------------
    fit = function(data, savespace = TRUE) {  
      assert_that(is.DatKeepClass(data))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) {
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$fit(data = data, savespace = savespace)
        }
        # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each fitRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        fitRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$fit(data = data, savespace = savespace)
        }
        # copy the fits one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.fit(fitRes[[k_i]])
        }
      }
      invisible(self)
    },  

    # take fitted GenericModel class object as an input and copy all the model fits down the line
    copy.fit = function(Generic.Model) {
      assert_that("GenericModel" %in% class(Generic.Model))
      private$PsAsW.models <- Generic.Model$getPsAsW.models()
      self$is.fitted <- TRUE
      invisible(self)
    },      
	      
    # P(A=1|W=w): uses private$m.fit to generate predictions
    predict = function(newdata, savespace = TRUE) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DatKeepClass(newdata))
      # serial loop over all regressions in PsAsW.models:
      if (!self$parfit_allowed) { 
        for (k_i in seq_along(private$PsAsW.models)) {
	  private$PsAsW.models[[k_i]]$predict(newdata = newdata, savespace = savespace)
	}
      # parallel loop over all regressions in PsAsW.models:
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = FALSE)
        # NOTE: Each predRes[[k_i]] will contain a copy of every single R6 object that was passed by reference ->
        # *** the size of fitRes is 100x the size of private$PsAsW.models ***
        predRes <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predict(newdata = newdata, savespace = savespace)
        }
        # copy the predictions one by one from BinOutModels above into private field for BinOutModels
        for (k_i in seq_along(private$PsAsW.models)) {
          private$PsAsW.models[[k_i]]$copy.predict(predRes[[k_i]])
        }
      }
      invisible(self)
    },

    # WARNING: This method cannot be chained together with other methods (s.a, class$predictAeqa()$fun())
    # Uses daughter objects (stored from prev call to fit()) to get predictions for P(A=obsdat.A|W=w)
    # Invisibly returns the joint probability P(A=a|W=w), also saves it as a private field "cumprodAeqa"
    # P(A=a|W=w) - calculating the likelihood for obsdat.A[i] (n vector of a's):
    predictAeqa = function(newdata, savespace = TRUE, wipeProb = TRUE, ...) {
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DatKeepClass(newdata))
      n <- newdata$nobs
      if (!self$parfit_allowed) {
	cumprodAeqa <- rep.int(1L, n)
        # loop over all regressions in PsAsW.models:
	for (k_i in seq_along(private$PsAsW.models)) {
	  cumprodAeqa <- cumprodAeqa * private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, savespace = savespace, wipeProb = wipeProb, ...)
	}
      } else if (self$parfit_allowed) {
        val <- checkpkgs(pkgs=c("foreach", "doParallel", "matrixStats"))
        mcoptions <- list(preschedule = TRUE)
        probAeqa_list <- foreach::foreach(k_i = seq_along(private$PsAsW.models), .options.multicore = mcoptions) %dopar% {
          private$PsAsW.models[[k_i]]$predictAeqa(newdata = newdata, savespace = savespace, wipeProb = wipeProb, ...)
        }
        # cbind_t <- system.time(
        probAeqa_mat <- do.call('cbind', probAeqa_list)
        # rowProds_t <- system.time(
        cumprodAeqa <- matrixStats::rowProds(probAeqa_mat)
      }
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  
  active = list(
    # recursively call all saved daughter model fits and wipe out any traces of saved data
    wipe.alldat = function(wipeProb = TRUE) {
      for (k_i in seq_along(private$PsAsW.models)) {
	private$PsAsW.models[[k_i]]$wipe.alldat <- wipeProb
      }
      return(self)
    }
  ),

  private = list(
    deep_clone = function(name, value) {
      # if value is is an environment, quick way to copy:
      # list2env(as.list.environment(value, all.names = TRUE), parent = emptyenv())
      # if a list of R6 objects, make a deep copy of each:
      if (name == "PsAsW.models") {
        lapply(value, function(PsAsW.model) PsAsW.model$clone(deep=TRUE))
        # to check the value is an R6 object:
      } else if (inherits(value, "R6")) {
        value$clone(deep=TRUE)
      } else {
        value  # For all other fields, just return the value  
      }
    },
    PsAsW.models = list(),
    fitted.pbins = list(),
    cumprodAeqa = NULL
  )
)
	      

# ---------------------------------------- def_regs_subset ----------------------------------
# Purpose: Define subset evaluation for new bins:
# same code in ContinModel$new and CategorModel$new replaced with outside function:
# -------------------------------------------------------------------------------------------
def_regs_subset <- function(self) {
  bin_regs <- self$reg$clone() # instead of defining new RegressionClass now cloning parent reg object and then ADDING new SETTINGS
  bin_regs$reg_hazard <- TRUE # don't add degenerate bins as predictors in each binary regression
  
  if (!self$reg$pool_cont) {
    bin_regs$outvar.class <- as.list(rep_len(gvars$sVartypes$bin, self$nbins))
    bin_regs$outvar <- self$bin_nms 
    bin_regs$predvars <- self$reg$predvars
    bin_regs$subset_vars <- lapply(self$bin_nms, function(var) { c(var, self$reg$subset_vars)})
    names(bin_regs$subset_vars) <- names(bin_regs$outvar.class) <- self$bin_nms
    
  } else {  # Same but when pooling across bin indicators:
    bin_regs$outvar.class <- gvars$sVartypes$bin
    bin_regs$outvar <- self$outvar
    bin_regs$outvars_to_pool <- self$bin_nms
    if (gvars$verbose)  {
      print("pooled bin_regs$outvar: "); print(bin_regs$outvar)
      print("bin_regs$outvars_to_pool: "); print(bin_regs$outvars_to_pool)
      print("bin_regs$subset_vars: "); print(bin_regs$subset_vars)
    }
  }
  bin_regs$resetS3class()
  return(bin_regs)
}
	      
	      
## -------------------------------------------------------------------------------------------
#' R6 class for modeling (fitting and predicting) joint probability for a univariate continuous outcome A[j]
#'
#'  \code{ContinModel} inherits from \code{\link{GenericModel}} class, defining and modeling a conditional density \code{P(A[j]|W,E,...)}
#'  where \code{A[j]} is univariate and continuous. By calling \code{self$new()}, \code{A[j]} will be discretized into \code{nbins} bins 
#'  via one of the 3 bin cutoff approaches (See Details for \code{\link{tmleCommunity}}). By calling \code{self$fit()}, it fits hazard 
#'  regressoin \code{Bin_A[j][t] ~ W + E} on \code{data} (a \code{\link{DatKeepClass}} class), which is the hazard probaility of the 
#'  the observation of A[j] belongs to bin \code{Bin_A[j][t]}, given covariates \code{(W, E)} and that observation doesn't belong to any
#'  precedent bins \code{Bin_A[j][1]}, \code{Bin_A[j][2]}, ..., \code{Bin_A[j][t-1]}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{nbins}} - 
#' \item{\code{bin_nms}} - Character vector of column names of bin indicators.
#' \item{\code{intrvls}} - 
#' \item{\code{intrvls.width}} - 
#' \item{\code{bin_weights}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, DataStorageClass.g0, DataStorageClass.gstar, ...)}}{Instantiate an new instance of \code{ContinModel} for a univariate continuous outcome A[j]}
#'   \item{\code{fit(data, savespace = TRUE)}}{...}
#'   \item{\code{predict(newdata, savespace = TRUE)}}{...}
#'   \item{\code{predictAeqa(newdata, savespace = TRUE, wipeProb = TRUE)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
ContinModel <- R6Class(classname = "ContinModel",
  inherit = GenericModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),        # The name of the continous outcome var (A[j])
    nbins = NULL,                # Actual nbins used, for cont. outvar
    bin_nms = character(),       # Column names for bin indicators
    intrvls = numeric(),         # Vector of numeric cutoffs defining the bins or a named list of numeric intervals (for length(self$outvar) > 1)
    intrvls.width = NULL,        # Named vector of bin-widths (bw_j : j=1,...,M) for each each bin in self$intrvls
    bin_weights = NULL,
    # Define settings for fitting contin A and then call $new for super class (GenericModel)
    initialize = function(reg, DatKeepClass.g0, DatKeepClass.gstar, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
 
      assert_that(is.DatKeepClass(DatKeepClass.g0))
      self$intrvls <- DatKeepClass.g0$detect.sVar.intrvls(reg$outvar,
                                                          nbins = self$reg$nbins,
                                                          bin_bymass = (getopt("bin.method") %in% "equal.mass"),
                                                          bin_bydhist = (getopt("bin.method") %in% "dhist"),
                                                          max_nperbin = as.integer(getopt("maxNperBin")))
      if (!missing(DatKeepClass.gstar)) {
        assert_that(is.DatKeepClass(DatKeepClass.gstar))
        gstar.intrvls <- DatKeepClass.gstar$detect.sVar.intrvls(reg$outvar,
                                                                nbins = self$reg$nbins,
                                                                bin_bymass = (getopt("bin.method") %in% "equal.mass"),
                                                                bin_bydhist = (getopt("bin.method") %in% "dhist"),
                                                                max_nperbin = as.integer(getopt("maxNperBin")))
        self$intrvls <- unique(sort(union(self$intrvls, gstar.intrvls)))
      }
        
      # Define the number of bins (no. of binary regressions to run),
      self$nbins <- self$reg$nbins <- length(self$intrvls) - 1
      # new outvar var names (bin names); all predvars remain unchanged;                                                    
      self$bin_nms <- DatKeepClass.g0$bin.nms.sVar(reg$outvar, self$nbins)
      # Save bin widths in reg class (naming the vector entries by bin names):
      self$intrvls.width <- diff(self$intrvls)
      self$intrvls.width[self$intrvls.width <= gvars$tolerr] <- 1
      self$reg$intrvls.width <- self$intrvls.width						  
      names(self$reg$intrvls.width) <- names(self$intrvls.width) <- self$bin_nms
      if (gvars$verbose)  {
        print("ContinModel outcome: " %+% self$outvar)
        print("ContinModel nbins: " %+% self$nbins)
      }
      bin_regs <- def_regs_subset(self = self)                                                          
      super$initialize(reg = bin_regs, ...)
    },

    # Transforms data for continous outcome to discretized bins A[j] -> BinsA[1], ..., BinsA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in A - names have changed though)    
    fit = function(data, savespace = TRUE) {
      assert_that(is.DatKeepClass(data))
      # Binirizes & saves binned matrix inside DatNet.sWsA
      data$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$nbins, bin.nms = self$bin_nms)
      if (gvars$verbose) {
        print("performing fitting for continuous outcome: " %+% self$outvar)
        print("freq counts by bin for continuous outcome: "); print(table(data$ord.sVar))
        print("binned dataset: "); print(head(cbind(data$ord.sVar, data$dat.bin.sVar), 5))
      }
      super$fit(data, savespace = savespace) # call the parent class fit method
      if (gvars$verbose) message("fit for outcome " %+% self$outvar %+% " succeeded...")
      if (savespace) data$emptydat.bin.sVar # wiping out binirized mat in data DatKeepClass object...
      if (savespace) self$wipe.alldat # wiping out all data traces in ContinModel...
      invisible(self)
    },
    # P(A=1|W=w): uses private$m.fit to generate predictions
    predict = function(newdata, savespace = TRUE) {
      if (gvars$verbose) print("performing prediction for continuous outcome: " %+% self$outvar)
      if (missing(newdata)) stop("must provide newdata")
      assert_that(is.DatKeepClass(newdata))      
      # mat_bin doesn't need to be saved (even though its invisibly returned); mat_bin is automatically saved in DatKeepClass - a potentially dangerous side-effect!!!
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$nbins, bin.nms = self$bin_nms)
      super$predict(newdata, savespace = savespace)
      if (savespace) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DatKeepClass object...
      invisible(self)
    },
    # Convert contin. A vector into matrix of binary cols, then call parent class method: super$predictAeqa()
    # Invisibly return cumm. prob P(A=a|W=w)
    predictAeqa = function(newdata, savespace = TRUE, wipeProb = TRUE) { # P(A=a|W=w) - calculating the likelihood for obsdat.A[i] (n vector of a's)
      if (gvars$verbose) print("performing prediction for continuous outcome: " %+% self$outvar)
      assert_that(is.DatKeepClass(newdata))
      newdata$binirize.sVar(name.sVar = self$outvar, intervals = self$intrvls, nbins = self$nbins, bin.nms = self$bin_nms)
      bws <- newdata$get.sVar.bw(name.sVar = self$outvar, intervals = self$intrvls)
      self$bin_weights <- (1 / bws) # weight based on 1 / (sVar bin widths)
      # Option 1: ADJUST FINAL PROB by bw.j TO OBTAIN density at a point f(sa|sw) = P(sA=sa|sW=sw):
      cumprodAeqa <- super$predictAeqa(newdata = newdata, savespace = savespace, wipeProb = wipeProb) * self$bin_weights
      # Alternative 2: ALso integrate the difference of sA value and its left most bin cutoff: x - b_{j-1} and pass it
      # This is done so that we can integrate the constant hazard all the way to the value of x:
        # * (1 - bw.j.sA_diff*(1/self$bin_weights)*probA1) (discrete)
        # * exp(-bw.j.sA_diff*(1/self$bin_weights)*probA1) (continuous)
      # bw.j.sA_diff <- newdata$get.sVar.bwdiff(name.sVar = self$outvar, intervals = self$intrvls)
      # cumprodAeqa <- super$predictAeqa(newdata = newdata, bw.j.sA_diff = bw.j.sA_diff) * self$bin_weights
      if (savespace) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      if (savespace) self$bin_weights <- NULL # wiping out self$bin_weights...
      if (savespace) self$wipe.alldat <- wipeProb # wiping out all data traces in ContinModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  active = list(
    cats = function() {seq_len(self$nbins)}
  )
)


## ---------------------------------------------------------------------
#' R6 class for modeling (fitting and predicting) joint probability for a univariate categorical outcome A[j]
#'
#'  \code{CategorModel} inherits from \code{\link{GenericModel}} class, defining and modeling a conditional density \code{P(A[j]|W,E...)}
#'  where \code{A[j]} is univariate and categorical. By calling \code{self$new()}, \code{A[j]} will be redefined into number of bins 
#'  \code{length(levels)} (i.e., number of unique categories in \code{A[j]}). By calling \code{self$fit()}, it fits hazard regressoin
#'  \code{Bin_A[j][t] ~ W + E} on \code{data} (a \code{\link{DatKeepClass}} class), which is the hazard probaility of the observation 
#'  of A[j] belongs to bin \code{Bin_A[j][t]}, given covariates \code{(W, E)} and that observation doesn't belong to any precedent bins 
#'  \code{Bin_A[j][1]}, \code{Bin_A[j][2]}, ..., \code{Bin_A[j][t-1]}.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#' \item{\code{reg}} - .
#' \item{\code{outvar}} - .
#' \item{\code{levels}} - Numeric vector of all unique categories in outcome outvar. 
#' \item{\code{nbins}} - .
#' \item{\code{bin_nms}} - .
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(reg, DatKeepClass.g0, ...)}}{Instantiate an new instance of \code{CategorModel} for a univariate categorical outcome A[j]}
#'   \item{\code{fit(data)}}{...}
#'   \item{\code{predict(newdata)}}{...}
#'   \item{\code{predictAeqa(newdata)}}{...}
#' }
#' @section Active Bindings:
#' \describe{
#'   \item{\code{cats}}{...}
#' }
#' @export
CategorModel <- R6Class(classname = "CategorModel",
  inherit = GenericModel,
  portable = TRUE,
  class = TRUE,
  public = list(
    reg = NULL,
    outvar = character(),     # the name of the categorical outcome var (A[j])
    levels = numeric(),       # all unique values for A[j] sorted in increasing order
    nbins = integer(),
    bin_nms = character(),
    # Define settings for fitting cat sA and then call $new for super class (SummariesModel)
    initialize = function(reg, DatKeepClass.g0, ...) {
      self$reg <- reg
      self$outvar <- reg$outvar
      # Define the number of bins (no. of binary regressions to run) based on number of unique levels for categorical sVar:
      # all predvars remain unchanged
      assert_that(is.DatKeepClass(DatKeepClass.g0))
      self$levels <- DatKeepClass.g0$detect.cat.sVar.levels(reg$outvar)

      self$nbins <- self$reg$nbins <- length(self$levels)
      self$bin_nms <- DatKeepClass.g0$bin.nms.sVar(reg$outvar, self$nbins)
      if (gvars$verbose)  {
        print("CategorSummaryModel outcome: "%+%self$outvar)
        print("CategorSummaryModel levels: "); print(self$levels)
      }
      bin_regs <- def_regs_subset(self = self)
      super$initialize(reg = bin_regs, ...)
    },
    # Transforms data for categorical outcome to bin indicators A[j] -> BinA[1], ..., BinA[M] and calls $super$fit on that transformed data
    # Gets passed redefined subsets that exclude degenerate Bins (prev subset is defined for names in A - names have changed though)
    fit = function(data, savespace = TRUE) {
      assert_that(is.DatKeepClass(data))
      # Binirizes & saves binned matrix inside DatKeepClass for categorical sVar
      data$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      if (gvars$verbose) {
        print("performing fitting for categorical outcome: " %+% self$outvar)
        print("freq counts by bin for categorical outcome: "); print(table(data$get.sVar(self$outvar)))
        print("binned dataset: "); print(head(cbind(sA = data$get.sVar(self$outvar), data$dat.bin.sVar), 5))
      }
      super$fit(data, savespace = savespace) # call the parent class fit method
      if (gvars$verbose) message("fit for " %+% self$outvar %+% " var succeeded...")
      if (savespace) data$emptydat.bin.sVar # wiping out binirized mat in data DatKeepClass object...
      # self$wipe.alldat # wiping out all data traces in CategorModel...
      invisible(self)
    },

    # P(A=1|W=w): uses private$m.fit to generate predictions
    predict = function(newdata, savespace = TRUE) {
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      if (missing(newdata)) stop("must provide newdata")  
      assert_that(is.DatKeepClass(newdata))
      newdata$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      super$predict(newdata, savespace = savespace)
      if (savespace) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata DatKeepClass object...
      invisible(self)
    },

    # Invisibly return cumm. prob P(A=a|W=w)
    # P(A=a|W=w) - calculating the likelihood for obsdat.A[i] (n vector of a's):
    predictAeqa = function(newdata, savespace = TRUE, wipeProb = TRUE) {
      if (gvars$verbose) print("performing prediction for categorical outcome: " %+% self$outvar)
      assert_that(is.DatKeepClass(newdata))
      newdata$binirize.cat.sVar(name.sVar = self$outvar, levels = self$levels)
      cumprodAeqa <- super$predictAeqa(newdata = newdata, savespace = savespace, wipeProb = wipeProb)
      if (savespace) newdata$emptydat.bin.sVar # wiping out binirized mat in newdata object...
      # self$wipe.alldat # wiping out all data traces in CategorModel...
      private$cumprodAeqa <- cumprodAeqa
      return(cumprodAeqa)
    }
  ),
  active = list(
    cats = function() {seq_len(self$nbins)}
  )
)
