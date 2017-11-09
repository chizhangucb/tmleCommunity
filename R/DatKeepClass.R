#----------------------------------------------------------
# General utilities / Global Vars
#----------------------------------------------------------
is.integerish <- function (x) { is.integer(x) || (is.numeric(x) && all(x == as.integer(x))) }

iqr <- function(x) { return(diff(quantile(x,c(.25, .75),na.rm=T))) }  # interquartile range

# Discarding missing outcomes (by default).
#------------------------------------
#' Panel Dataset Transformation, Using Individual (and Time) Indexes
#'
#' \code{panelData_Trans} provides a wide variety of ways of data transformation for panel datasets, such as fixid 
#'  effect and pooling model. It allows users to only apply transformation on regressors of interests, instead of
#'  on the entire dataset. See details in \url{https://cran.r-project.org/web/packages/plm/plm.pdf}.
#   and \url{https://github.com/cran/plm/blob/master/R/pFormula.R} & \url{https://github.com/cran/plm/blob/master/R/plm.R}
#' @param data A data frame (will be automatically transferred to panel data frame) or a panel data frame
#' @param yvar Column name in \code{data} of outcome variable (Currently only support univariate).    
#' @param xvar Column names in \code{data} of explanatory variables (Including \eqn{(A, W, E)}).
#' @param effect The effects introduced in the model, one of "individual", "time", "twoways" and "nested". Default to "individual".
#' @param model Model of estimation, one of "pooling" (pooled OLS), "within" (fixed effect), "between" (group mean), 
#'  "random"(random effect), "fd" (first differences) and "ht" (Hausman-Taylor estimator). Default to "within". Notice that when 
#'  \code{model} = "random", "swar" is chosen as the method of estimation for the variance components
#' @param index A vector of two character strings which contains the names of the individual and of the time indices, sequentially. 
#'  If only individual index is given, treat each observation within a unit as a time point.
#'  If no index is given, the first two columns will be automatically treated as individual and time indices, sequentially.
#' @param transY Logical. If \code{TRUE}, indicate the outcome variable \code{yvar} will also be tranformed. Default to \code{TRUE}.
#' @return \item \code{newdata} Transformed panel data
#' @example tests/examples/8_panelData_Trans_examples.R
#' @export
panelData_Trans <- function(data, yvar, xvar, effect = "individual", model = "within", index = NULL, transY = TRUE) {
  formula <- as.formula(yvar %+% " ~ " %+% paste(xvar, collapse=" + "))
  # Check whether data is a pdata.frame and if not create it
  orig_rownames <- row.names(data)
  
  # Create a data.frame with an index attribute that describes its individual and time dimensions
  if (!inherits(data, "pdata.frame")) data <- plm::pdata.frame(data, index)
  # Create a pFormula object if necessary
  if (!inherits(formula, "pFormula")) formula <- plm::pFormula(formula)
  
  # eval the model.frame (Nothing to do with the trnsformation)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf$formula <- formula
  mf$data <- data
  # eval in parent.frame() doesn't work
  #  data <- eval(mf, sys.frame(which = nframe))
  data <- eval(mf, parent.frame())  # It drops the columns that are not indicated in the formula
  
  # preserve original row.names for data [also fancy rownames]; so functions like
  # pmodel.response(), model.frame(), model.matrix(), residuals() return the original row.names
  # Reason: eval(mf, parent.frame()) returns row.names as character vector containing 
  # "row_number" with incomplete observations dropped
  row.names(data) <- orig_rownames[as.numeric(row.names(data))]
  # return the model.frame or estimate the model
  if (is.na(model)){
    attr(data, "formula") <- formula
    return(data)
  }
  
  args <- list(model = model, effect = effect)
  
  # extract the model.matrix and the model.response (Transformation related)
  if (model == "random") {
    theta <- plm::ercomp(formula, data, index = index, effect = effect, method = "swar")$theta
  } else {
    theta <- NULL
  }
  X <- plm:::model.matrix.pFormula(object = formula, data, model = model, effect = effect, rhs = 1, theta = theta)
  if (ncol(X) == 0) stop("empty model")
  if (transY) {
    Y <- plm:::pmodel.response(object = formula, data, model = model, effect = effect, theta = theta)
  } else {
    Y <- data[, yvar]
  }
  
  newdata <- as.data.frame(cbind(Y, X))
  colnames(newdata) <- c(yvar, colnames(X))
  return(newdata)
}

## Original method code from https://github.com/hadley/densityvis/blob/master/R/breaks-dhist.r
## Modify code from https://github.com/cran/tmlenet/blob/master/R/dhist.r
dhist <- function(x, a = 5*iqr(x), nbins=nclass.Sturges(x), rx = range(x, na.rm=TRUE), eps=.15, xlab = "x", plot = TRUE, lab.spikes=TRUE) {
# Arguments:
#    x          - the data
#    a          - the scaling factor, default is 5 * IQR
#    nbins      - the number of bins, default is assigned by the Stuges method (basing bin sizes on the range of the data)
#    rx         - the range used for the left of the left-most bin to the right of the right-most bin
#    xlab       - label for the x axis
#    plot       - Logical, if TRUE produces the plot, if FALSE returns the heights, breaks and counts
#    lab.spikes - Logical, if TRUE labels the % of data in the spikes
if (is.character(nbins)) {
    nbins <- switch(casefold(nbins),
                    sturges = nclass.Sturges(x),
                    fd = nclass.FD(x),
                    scott = nclass.scott(x),
                    stop("Nclass method not recognized"))
  } else if (is.function(nbins)) {
    nbins <- nbins(x)
  }
  x <- sort(x[!is.na(x)])
  
  if (a == 0) { a <- diff(range(x))/100000000 }
  if (a != 0 & a != Inf) {
    n <- length(x)
    h <- (rx[2] + a - rx[1]) / nbins
    ybr <- rx[1] + h * (0:nbins)
    yupper <- x + (a * (1:n))/n
    # upper and lower corners in the ecdf
    ylower <- yupper - a/n
    
    cmtx <- cbind(cut(yupper, breaks = ybr), cut(yupper, breaks = ybr, left.include = T), 
                  cut(ylower, breaks = ybr), cut(ylower, breaks = ybr, left.include = T))
    cmtx[1, 3] <- cmtx[1, 4] <- 1
    # to replace NAs when default r is used
    cmtx[n, 1] <- cmtx[n, 2] <- nbins
    # checksum <- apply(cmtx, 1, sum) %% 4
    checksum <- (cmtx[, 1] + cmtx[, 2] + cmtx[, 3] + cmtx[, 4]) %% 4
    # will be 2 for obs. that straddle two bins
    straddlers <- (1:n)[checksum == 2]
    # to allow for zero counts
    if(length(straddlers) > 0) {
      counts <- table(c(1:nbins, cmtx[- straddlers, 1])) 
    } else {
      counts <- table(c(1:nbins, cmtx[, 1]))
    }
    counts <- counts - 1
    
    if(length(straddlers) > 0) {
      for(i in straddlers) {
        binno <- cmtx[i, 1]
        theta <- ((yupper[i] - ybr[binno]) * n)/a
        counts[binno - 1] <- counts[binno - 1] + (1 - theta)
        counts[binno] <- counts[binno] + theta
      }
    }
    xbr <- ybr
    xbr[-1] <- ybr[-1] - (a * cumsum(counts))/n
    spike <- eps * diff(rx) / nbins
    flag.vec <- c(diff(xbr) < spike, FALSE)
    if ( sum(abs(diff(xbr)) <= spike) > 1 ) {
      xbr.new <- xbr
      counts.new <- counts
      diff.xbr <- abs(diff(xbr))
      amt.spike <- diff.xbr[length(diff.xbr)]
      for (i in rev(2:length(diff.xbr))) {
        if ((diff.xbr[i-1] <= spike) & (diff.xbr[i] <= spike) & (!is.na(diff.xbr[i]))) {
          amt.spike <- amt.spike + diff.xbr[i-1]
          counts.new[i-1] <- counts.new[i-1] + counts.new[i]
          xbr.new[i] <- NA
          counts.new[i] <- NA
          flag.vec[i-1] <- T
        } else {
          amt.spike <- diff.xbr[i-1]
        }
      }
      flag.vec <- flag.vec[!is.na(xbr.new)]
      flag.vec <- flag.vec[-length(flag.vec)]
      counts <- counts.new[!is.na(counts.new)]
      xbr <- xbr.new[!is.na(xbr.new)]
    } else {
      flag.vec<-flag.vec[-length(flag.vec)]
    }
    widths <- abs(diff(xbr))
    # N.B. argument "widths" in barplot must be xbr
    heights <- counts / widths
  }
  bin.size <- length(x) / nbins
  cut.pt <- unique(c(min(x) - abs(min(x))/1000, 
                     approx(seq(length(x)), x, 
                            (1:(nbins - 1)) * bin.size, rule = 2)$y, 
                     max(x)))
  aa <- hist(x, breaks = cut.pt, plot = FALSE)
  # aa <- hist(x, breaks = cut.pt, plot = FALSE, probability = TRUE)
  if(a == Inf) {
    heights <- aa$counts
    xbr <- aa$breaks
  }
  amt.height <- 3
  q75 <- quantile(heights, .75)
  if (sum(flag.vec) != 0) {
    amt <- max(heights[!flag.vec])
    ylim.height <- amt * amt.height
    ind.h <- flag.vec&heights > ylim.height
    flag.vec[heights < (ylim.height * (amt.height - 1) / amt.height)] <- FALSE
    heights[ind.h] <- ylim.height
  }
  amt.txt <- 0
  end.y <- (-10000)
  if (plot) {
    barplot(heights, abs(diff(xbr)), space = 0, density = -1, 
            xlab = xlab, plot = TRUE, xaxt = "n",yaxt='n')
    at <- pretty(xbr)
    axis(1, at = at - xbr[1], labels = as.character(at))
    if (lab.spikes) { 
      if (sum(flag.vec) >= 1) {
        usr <- par('usr')
        for ( i in seq(length(xbr) - 1)) {
          if (!flag.vec[i]) {
            amt.txt <- 0
            if ((xbr[i] - xbr[1]) < end.y) amt.txt <- 1
          }
          else {
            amt.txt <- amt.txt + 1
            end.y <- xbr[i] - xbr[1] + 3 * par('cxy')[1]
          }
          if (flag.vec[i]) {
            txt< - paste(' ', format(round(counts[i] / sum(counts)*100)), '%', sep='')
            par(xpd = T)
            text(xbr[i+1] - xbr[1], ylim.height - par('cxy')[2] * (amt.txt - 1), txt, adj=0)
          }}
      } else {
        print('no spikes or more than one spike')
      }
    }
    invisible(list(heights = heights, xbr = xbr))
  } else {
    return(list(heights = heights, xbr = xbr, counts = counts))
  }
}

# Get the actual A^* sampled from user-supplied intervention f_gstar (function_name):
f.gen.A.star <- function(data, f.g_fun) {
  .f_g_wrapper <- function(data, f.g_fun, ...) {
    args0 <- list(data = data)
    args <- c(args0, ...)
    do.call(f.g_fun, args)
  }
  # test f.g_fun is a function, if not it must be a vector or a matrix/ data.frame
  if (!is.function(f.g_fun)) {
    if (is.matrix(f.g_fun) || is.data.frame(f.g_fun)) {
      newA <- as.data.frame(f.g_fun)
      if (nrow(newA)!=nrow(data) && nrow(newA)!=1L) 
        stop("f_gstar1/f_gstar2 must be either a function or a vector of length nrow(data) or 1  
             or a data frame/ matrix with nrow(data) or 1")
      if (nrow(newA)==1L) newA <- do.call(cbind, lapply(1:ncol(newA), function(i) rep_len(newA[1, i], nrow(data))))
    } else if (is.vector(f.g_fun)) {
      newA <- as.vector(f.g_fun)
      if (length(newA) != nrow(data) && length(newA) != 1L) {
        stop("f_gstar1/f_gstar2 must be either a function or a vector of length nrow(data) or 1 
             or a data frame/ matrix of nrow(data) or 1 row(s)")
      }
      if (length(newA)==1L) newA <- rep_len(newA, nrow(data))      
    }
  } else {
    if (!("data" %in% names(formals(f.g_fun)))) 
      stop("functions f_gstar1 / f_gstar2 must have a named argument 'data'")  
    newA <- .f_g_wrapper(data = data, f.g_fun = f.g_fun)
  }
  return(newA)
}

## ---------------------------------------------------------------------
# DETECTING VECTOR TYPES
# sVartypes <- list(bin = "binary", cat = "categor", cont = "contin")
## ---------------------------------------------------------------------
detect.col.types <- function(sVar_mat){
  assert_that(is.integerish(getopt("maxncats")) && getopt("maxncats") > 1)
  maxncats <- getopt("maxncats")
  sVartypes <- gvars$sVartypes
  as.list(apply(sVar_mat, 2, function(vec) {
    vec_nomiss <- vec[!gvars$misfun(vec)]
    nvals <- length(unique(vec_nomiss))
    if (nvals <= 2L) { sVartypes$bin
    } else if ((nvals <= maxncats) && (is.integerish(vec_nomiss))) { sVartypes$cat
    } else { sVartypes$cont }
  }))
}

## ---------------------------------------------------------------------
# Normalizing / Defining bin intervals / Converting contin. to ordinal / Converting ordinal to bin indicators
## ---------------------------------------------------------------------
normalize <- function(x) {
  if (abs(max(x) - min(x)) > gvars$tolerr) { # Normalize to 0-1 only when x is not constant
    return((x - min(x)) / (max(x) - min(x)))
  } else {  # What is the thing to do when x constant? Set to abs(x), abs(x)/x or 0???
    return(x)
  }
}

# Define bin cutt-offs for continuous x:
define.intervals <- function(x, nbins, bin_bymass, bin_bydhist, max_nperbin) {
  x <- x[!gvars$misfun(x)]  # remove missing vals
  nvals <- length(unique(x))
  if (is.na(nbins)) nbins <- as.integer(length(x) / max_nperbin)
  # if nbins is too high, for ordinal, set nbins to n unique obs and cancel quantile based interval defns
  if (nvals < nbins) {
    nbins <- nvals
    bin_bymass <- FALSE
  }
  if (abs(max(x) - min(x)) > gvars$tolerr) {  # when x is not constant
    if ((bin_bymass) & !is.null(max_nperbin)) {
      if ((length(x) / max_nperbin) > nbins) nbins <- as.integer(length(x) / max_nperbin)
    }
    intvec <- seq.int(from = min(x), to = max(x) + 1, length.out = (nbins + 1)) # interval type 1: bin x by equal length intervals of 0-1
  } else {  # when x is constant, force the smallest possible interval to be at least [0,1]
    intvec <- seq.int(from = min(0L, min(x)), to = max(1L, max(x)), length.out = (nbins + 1))
  }
  if (bin_bymass) {
    intvec <- quantile(x = x, probs = normalize(intvec)) # interval type 2: bin x by mass (quantiles of 0-1 intvec as probs)
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  } else if (bin_bydhist) {
    # stop("... binning continuous variable: dhist bin definitions are no longer available...")
    intvec <- dhist(x, plot = FALSE, nbins = nbins)$xbr
    intvec[1] <- intvec[1] - 0.01
    intvec[length(intvec)] <- intvec[length(intvec)] + 0.01
  }
  # adding -Inf & +Inf as leftmost & rightmost cutoff points to make sure all future data points end up in one of the intervals:
  intvec <- c(min(intvec)-1000L, intvec, max(intvec)+1000L)
  return(intvec)
}

# Turn any x into ordinal (1, 2, 3, ..., nbins) for a given interval cutoffs (length(intervals)=nbins+1)
make.ordinal <- function(x, intervals) findInterval(x = x, vec = intervals, rightmost.closed = TRUE)

# Make dummy indicators for ordinal x (sA[j]) Approach used: creates B_j that jumps to 1 only once
# and stays 1 (degenerate) excludes reference category (last)                                                       
make.bins_mtx_1 <- function(x.ordinal, nbins, bin.nms, levels = 1:nbins) {
  n <- length(x.ordinal)
  new.cats <- 1:nbins
  dummies_mat <- matrix(1L, nrow = n, ncol = length(new.cats))
  for(cat in new.cats[-length(new.cats)]) {
    subset_Bj0 <- x.ordinal > levels[cat]
    dummies_mat[subset_Bj0, cat] <- 0L
    subset_Bjmiss <- x.ordinal < levels[cat]
    dummies_mat[subset_Bjmiss, cat] <- gvars$misval
  }
  dummies_mat[, new.cats[length(new.cats)]] <- gvars$misval
  colnames(dummies_mat) <- bin.nms
  return(dummies_mat)
}


#-----------------------------------------------------------------------------
# Class Membership Tests
#-----------------------------------------------------------------------------
is.DatKeepClass <- function(DatKeepClass) "DatKeepClass" %in% class(DatKeepClass)


## ---------------------------------------------------------------------
#' R6 class for Storing, Managing, Subsetting and Manipulating the Input Data.
#'
#'  \code{DatKeepClass} allows user to access the input data. The processed covariates from sVar.object are stored as a matrix in
#'  (\code{private$.mat.sVar}). This class could subset, combine, normalize, discretize and binarize covariates in (A, W, E). 
#'  For disretization of continous and categorical variables, it can automatically detect / set covariates type (binary, categor, 
#'  contin), detect / set bin intervals, and construct bin indicators. Besides, it provides methods for generating new exposures under 
#'  user-specific arbitrary intervention \eqn{g^{*}} through \code{self$make.dat.sVar}, and allows user to replace missing values with                                                 
#'  user-specific \code{gvars$misXreplace} (Default to 0). Its pointers will be passed on to \code{GenericModel} functions: using in
#'  \code{$fit()}, \code{$predict()} and \code{$predictAeqa()}.
#'  
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#' @details
#' \itemize{
#'    \item{\code{norm.c.sVars}} - \code{flag} if \code{TRUE} normalize continous covariates.
#'    \item{\code{mat.bin.sVar}} - Matrix of the binary indicators created from discretization of continuous covariate \code{active.bin.sVar}.
#'    \item{\code{ord.sVar}} - Ordinal (categorical) transformation of a continous covariate \code{sVar}.
#'    \item{\code{obs.wts}} - Vectopr of observation (sampling) weights (of length \code{ndat.sVar}). If NULL, assumed to be all 1. 
#'    \item{\code{YnodeVals}} - Vector of outcome values (Ynode) in observed data
#'    \item{\code{det.Y}} - Logical vector, where \code{YnodeVals[det.Y==TRUE]} are deterministic and set to NA. 
#'    \item{\code{p}} - Number of Monte-Carlo simulations performed. 
#'    \item{\code{ndat.sVar}} - Number of observations in the observed data frame.
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(Odata, nodes, YnodeVals, det.Y, norm.c.sVars = FALSE, ..)}}{Instantiate an new instance of \code{DatKeepClass} that is 
#'     used for storing and manipulating the input data.}
#'   \item{\code{addYnode(YnodeVals, det.Y)}}{Add protected Y node to private field and set to NA all determinisitc Y values for public field YnodeVals.}
#'   \item{\code{addObsWeights(obs.wts)}}{Add observation weights to public field.}                                                        
#'   \item{\code{evalsubst(subset_vars, subset_exprs = NULL)}}{...}
#'   \item{\code{get.dat.sVar(rowsubset = TRUE, covars)}}{Subset covariate design matrix for \code{\link{BinaryOutModel}}.}
#'   \item{\code{get.outvar(rowsubset = TRUE, var)}}{Subset a vector of outcome variable for \code{BinaryOutModel}.}
#'   \item{\code{get.obsweights(rowsubset = TRUE)}}{Subset a vector of observation weights for \code{BinaryOutModel}.}
#'   \item{\code{def.types.sVar(type.sVar = NULL)}}{Define each variable' class in input data: bin, cat or cont.}
#'   \item{\code{set.sVar.type(name.sVar, new.type)}}{Assign a new class type to one variable that belongs to the input data.}
#'   \item{\code{get.sVar.type(name.sVar)}}{Return the class type of a variable.}
#'   \item{\code{is.sVar.cont(name.sVar)}}{Check if the variable is continuous.}
#'   \item{\code{is.sVar.cat(name.sVar)}}{Check if the variable is categorical.}
#'   \item{\code{is.sVar.bin(name.sVar)}}{Check if the variable is binary.}
#'   \item{\code{get.sVar(name.sVar)}}{Return a vector of the variable values.}
#'   \item{\code{set.sVar(name.sVar, new.sVarVal)}}{Assign a vector of new values to the specific variable.}
#'   \item{\code{bin.nms.sVar(name.sVar, nbins)}}{Define names of bin indicators for \code{sVar.}}
#   \item{\code{pooled.bin.nm.sVar(name.sVar)}}{Define a name of pooled bin indicators for \code{sVar}.}
#'   \item{\code{detect.sVar.intrvls(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin)}}{...}
#'   \item{\code{detect.cat.sVar.levels(name.sVar)}}{Detect the unique categories in categorical sVar, returning in increasing order.}
#'   \item{\code{get.sVar.bw(name.sVar, intervals)}}{Get the bin widths vector for the discretized cont \code{sVar}.}
#'   \item{\code{get.sVar.bwdiff(name.sVar, intervals)}}{Get the bin widths differences vector for the discretized continuous \code{sVar}.}
#'   \item{\code{binirize.sVar(name.sVar, ...)}}{Create a matrix of bin indicators for categorical/cont \code{sVar}.}
#'   \item{\code{norm.cont.sVars()}}{Normalize continuous \code{sVar} (Note that this process is memory-costly).}
#'   \item{\code{fixmiss_sVar()}}{Replace all missing (NA) values with a default integer (Default to 0).}
#'   \item{\code{make.dat.sVar(p = 1, f.g_fun = NULL, regform = NULL)}}{Generate new exposures under user-specific arbitrary intervention 
#'     \code{f.g_fun} and construct a data.frames that combines all covariates, replacing the old exposures with the new ones.}
#' }
#' @section Active Bindings:
#' \describe{
#'    \item{\code{names.sVar}}{Return variable names of the input data.}
#'    \item{\code{names.c.sVar}}{Return continuous variable names of the input data.}
#'    \item{\code{ncols.sVar}}{Return the number of columns of the input data.}
#'    \item{\code{nobs}}{Return the number of observations of the input data.}
#'    \item{\code{dat.sVar}}{Return a data frame object that stores the entire dataset (including all \code{sVar}.).}
#'    \item{\code{dat.bin.sVar}}{Return a stored matrix for bin indicators on currently binarized continous \code{sVar}.}
#'    \item{\code{active.bin.sVar}}
#'      {Return name(s) of active binarized continous sVar(s), changing when \code{fit} or \code{predict} is called.}
#'    \item{\code{emptydat.sVar}}{Wipe out \code{dat.sVar}.}
#'    \item{\code{emptydat.bin.sVar}}{Wipe out \code{dat.bin.sVar}.}
#'    \item{\code{noNA.Ynodevals}}{Return the observed Y without any missing values.}
#'    \item{\code{nodes}}{...}
#'    \item{\code{type.sVar}}{...}
#' }
#' @importFrom assertthat assert_that is.count
#' @seealso \code{\link{tmleCom_Options}}, \code{\link{tmleCommunity}}
#' @example tests/examples/2_DatKeepClass_examples.R                                                        
#' @export
DatKeepClass <- R6Class(classname = "DatKeepClass",
  portable = TRUE,
  class = TRUE,
  public = list(
    norm.c.sVars = FALSE,      # flag = TRUE if want to normalize continous covariates
    mat.bin.sVar = NULL,       # temp storage mat for bin indicators on currently binarized continous sVar (from self$active.bin.sVar)
    ord.sVar = NULL,           # Ordinal (cat) transform for continous sVar
    obs.wts = NULL,            # Observation (sampling) weights
    YnodeVals = NULL,          # Values of the binary outcome (Ynode) in observed data where det.Y = TRUE obs are set to NA
    det.Y = NULL,              # Logical vector, where YnodeVals[det.Y==TRUE] are deterministic (0 or 1)
    p = 1,                     # Number of Monte-Carlo simulations performed 
    ndat.sVar = NA_integer_,   # n of samples in the OBSERVED (original) data

    initialize = function(Odata, nodes, YnodeVals, det.Y, obs.wts, norm.c.sVars = FALSE, ...) {
      ######## CAUTIONS ########
      ## Not sure if want to use data.table all the time, so decide to keep as data frame first!
      assert_that(is.data.frame(Odata))  # | is.data.table(Odata))
      if (!is.list(nodes) || any(!(names(nodes) %in% c("Ynode", "Anodes", "WEnodes", "Crossnodes", "communityID"))) ) { 
        message("Don't recognize " %+% paste0(setdiff(names(nodes), c("Ynode", "Anodes", "WEnodes", "communityID", "Crossnodes")), collapse = " and "))
        stop("It should be a list & its names can only be one or more of Ynode, Anodes, WEnodes, communityID and Crossnodes.\n")
      } else { 
        self$nodes <- nodes 
      }
      if (!missing(YnodeVals)) self$addYnode(YnodeVals = YnodeVals, det.Y = det.Y)
      if (!missing(obs.wts)) self$addObsWeights(obs.wts = obs.wts)
        
      self$dat.sVar <- as.data.frame(Odata[, Reduce(c, nodes)])  # as.data.table(Odata)  # dat.sVar - makes a copy of the input data (shallow)
      self$ndat.sVar <- NROW(Odata)
      if (norm.c.sVars) {
        self$norm.c.sVars <- norm.c.sVars
        self$norm.cont.sVars()
      }  
      if (gvars$verbose) cat("...detecting the type of each input column...")
      self$def.types.sVar() # Define the type of each sVar[i]: bin, cat or cont

      invisible(self)
    },
        
    # -------------------------------------- addYnode -----------------------------------
    # Purpose: Add protected Y nodes to private field and set to NA all determinisitc Y values for public field YnodeVals
    # -----------------------------------------------------------------------------------
    addYnode = function(YnodeVals, det.Y) {
      if (missing(det.Y)) det.Y <- rep.int(FALSE, length(YnodeVals))
      self$noNA.Ynodevals <- YnodeVals  # Adding actual observed Y as protected (without NAs)
      self$YnodeVals <- YnodeVals
      self$YnodeVals[det.Y] <- NA       # Adding public YnodeVals & setting det.Y values to NA
      self$det.Y <- det.Y
    },
      
    # ------------------------------------ addObsWeights --------------------------------
    # Purpose: Add observation weights to public field  
    # -----------------------------------------------------------------------------------
    addObsWeights = function(obs.wts) {
      if (any(obs.wts < 0)) { stop("Observation weights have to be non-negative") }
      if (length(obs.wts) == 1) {
        self$obs.wts <- rep(obs.wts, self$ndat.sVar) 
      } else if (length(obs.wts) != self$ndat.sVar) { 
        stop("The length of observation weights should be the same as nrow(data)") 
      }
      self$obs.wts <- obs.wts
    },

    # ----------------------------------- evalsubst -------------------------------------
    # Purpose: Eval the expression (in the environment of the data.frame "data" + global constants "gvars")
    # -----------------------------------------------------------------------------------
    # Could also do evaluation in a special env with a custom subsetting fun '[' that will dynamically find the 
    # correct dataset that contains sVar.name (dat.sVar or dat.bin.sVar) and will return sVar vector
    evalsubst = function(subset_vars, subset_exprs = NULL) {
      if (!missing(subset_vars)) {
        assert_that(is.character(subset_vars))
        res <- rep.int(TRUE, self$nobs)
        for (subsetvar in subset_vars) {
          # (*) find the var of interest (in self$dat.sVar or self$dat.bin.sVar), give error if not found
          sVar.vec <- self$get.outvar(var = subsetvar)
          assert_that(!is.null(sVar.vec))
          # (*) reconstruct correct expression that tests for missing values
          res <- res & (!gvars$misfun(sVar.vec))
        }
        return(res)
      }
      if (!is.null(subset_exprs)) {
        if (is.logical(subset_exprs)) {
          return(subset_exprs)
        } else if (is.character(subset_exprs)) {
          ## ******************************************************
          ## data.table evaluation of the logical subset expression. May have 500K w/ 1000 bins: 4-5sec
          ## Note: This can be made faster by using keys in data.table on variables in eval(parse(text = subset_exprs))
          ## ******************************************************
          eval.env <- c(data.frame(self$dat.sVar), data.frame(self$dat.bin.sVar), as.list(gvars))
          res <- try(eval(subset_exprs, envir = eval.env, enclos = baseenv()))
          return(res)
        } else if (is.integer(subset_exprs)) {  
          ## The expression is a row index, hence should be returned to logical
          res <- rep.int(FALSE, self$nobs)
          res[subset_exprs] <- TRUE
          return(res)
        }
      }
    },

    # ----------------------------------- get.dat.sVar ------------------------------------
    # Purpose: subsetting/returning covariate design mat for BinaryOutModel Class
    # -------------------------------------------------------------------------------------
    get.dat.sVar = function(rowsubset = TRUE, covars) {
      if (!missing(covars)) {
        if (length(unique(colnames(self$dat.sVar))) < length(colnames(self$dat.sVar))) {
          warning("repeating column names in the final data set; please check for duplicate covariate / node names")
        }
        # columns to select from main design matrix (in the same order as listed in covars):
        sel.sVar <- intersect(covars, colnames(self$dat.sVar))
        if (is.matrix(self$dat.sVar) | is.data.frame(self$dat.sVar) ) {
          dfsel <- self$dat.sVar[rowsubset, sel.sVar, drop = FALSE] # data stored as matrix or data frame
        } else {
          stop("self$dat.sVar is of unrecognized class: " %+% class(self$dat.sVar))
        } # else if (is.data.table(self$dat.sVar)) {
          # dfsel <- self$dat.sVar[rowsubset, sel.sVar, drop = FALSE, with = FALSE] # data stored as data.table
          # } 
        # columns to select from binned continuous/cat var matrix (if it was previously constructed):
        if (!is.null(self$dat.bin.sVar)) {
          sel.binsVar <- intersect(covars, colnames(self$dat.bin.sVar))
        } else {
          sel.binsVar <- NULL
        }
        if (length(sel.binsVar)>0) dfsel <- cbind(dfsel, self$dat.bin.sVar[rowsubset, sel.binsVar, drop = FALSE])
         
        found_vars <- covars %in% colnames(dfsel)  # checks for non-existance of a particular var in covars
        if (!all(found_vars)) stop("some covariates can't be found: "%+% paste(covars[!found_vars], collapse=","))
        return(dfsel)
      } else {
        return(self$dat.sVar[rowsubset, , drop = FALSE])
      }
    },

    # ------------------------------------- get.outvar ------------------------------------
    # Purpose: subsetting/returning a vector of outcome variable for BinaryOutModel Class
    # -------------------------------------------------------------------------------------
    get.outvar = function(rowsubset = TRUE, var) {
      if (length(self$nodes) < 1) stop("DatKeepClass$nodes list is empty!")
      if (var %in% self$names.sVar) {
        self$dat.sVar[rowsubset, var]
        # if (is.data.table(self$dat.sVar)) {
        #   self$dat.sVar[rowsubset, var, with = FALSE]
        # } else {
        #   self$dat.sVar[rowsubset, var, drop = FALSE]
        # }
      } else if (var %in% colnames(self$dat.bin.sVar)) {
        self$dat.bin.sVar[rowsubset, var]
      } else if ((var %in% self$nodes$Ynode) && !is.null(self$YnodeVals)) {
        self$YnodeVals[rowsubset]
      } else {
        stop("requested variable " %+% var %+% " does not exist in DatKeepClass!")
      }
    },
      
    # ----------------------------------- get.obsweights ----------------------------------
    # Purpose: subsetting/returning a vector of observation weights for BinaryOutModel Class
    # -------------------------------------------------------------------------------------
    get.obsweights = function(rowsubset = TRUE) {
      if (!is.null(self$obs.wts)) {
        self$obs.wts[rowsubset]
      } else {
        rep(1, self$ndat.sVar)[rowsubset]
      }
    },

    # --------------------------------------------------------------------------------------
    # Methods for sVar types - Define the type (class) of each variable (sVar) in input data: bin, cat or cont
    # 5 ancillary functions: def.types.sVar, set.sVar.type, get.sVar.type, is.sVar.cont, is.sVar.cat
    # --------------------------------------------------------------------------------------
    # type.sVar acts as a flag: only detect types when !is.null(type.sVar), 
    # otherwise can pass type.sVar = list(sVar = NA, ...) or a value type.sVar = NA/gvars$sVartypes$bin/etc
    def.types.sVar = function(type.sVar = NULL) {
      if (is.null(type.sVar)) {
        private$.type.sVar <- detect.col.types(self$dat.sVar)
      } else {
        n.sVar <- length(self$names.sVar)
        len <- length(type.sVar)
        assert_that((len == n.sVar) || (len == 1L))
        if (len == n.sVar) { # set types for each variable
          assert_that(is.list(type.sVar))
          assert_that(all(names(type.sVar) %in% self$names.sVar))
        } else { # set one type for all vars
          assert_that(is.string(type.sVar))
          type.sVar <- as.list(rep(type.sVar, n.sVar))
          names(type.sVar) <- self$names.sVar
        }
        private$.type.sVar <- type.sVar
      }
      invisible(self)
    },

    set.sVar.type = function(name.sVar, new.type) { private$.type.sVar[[name.sVar]] <- new.type },
    get.sVar.type = function(name.sVar) { 
      if (missing(name.sVar)) { private$.type.sVar } else { private$.type.sVar[[name.sVar]] } 
    },    
    is.sVar.cont = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cont },    
    is.sVar.cat = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$cat },
    is.sVar.bin = function(name.sVar) { self$get.sVar.type(name.sVar) %in% gvars$sVartypes$bin },
      
    # --------------------------------------------------------------------------------------
    # Methods for directly handling one continous/categorical sVar in self$mat.sVar;
    # --------------------------------------------------------------------------------------
    get.sVar = function(name.sVar) { # self$dat.sVar[, name.sVar, with=FALSE]
      x <- self$dat.sVar[, name.sVar]
      if (is.list(x) || is.data.frame(x)) x <- x[[1]]
      return(x)
    },
    set.sVar = function(name.sVar, new.sVarVal) {
      assert_that(length(new.sVarVal)==self$nobs | length(new.sVarVal)==1)
      assert_that(name.sVar %in% colnames(self$dat.sVar))
      self$dat.sVar[, name.sVar] <- new.sVarVal  
      # self$dat.sVar[, (name.sVar) := new.sVarVal] if data.table
      invisible(self)
    },  
  
    # --------------------------------------------------------------------------------------
    # Binning methods for categorical/continuous sVar
    # 7 ancillary functions: bin.nms.sVar, pooled.bin.nm.sVar, detect.sVar.intrvls, 
    #                        detect.cat.sVar.levels, get.sVar.bw, get.sVar.bwdiff
    # --------------------------------------------------------------------------------------
    # Need to find a way to over-ride nbins for categorical vars (allowing it to be set to more than gvars$maxncats)!
    # Return names of bin indicators for sVar:
    bin.nms.sVar = function(name.sVar, nbins) { name.sVar %+% "_" %+% "B." %+% (1L:nbins) },
    pooled.bin.nm.sVar = function(name.sVar) { name.sVar %+% "_allB.j" },
    detect.sVar.intrvls = function(name.sVar, nbins, bin_bymass, bin_bydhist, max_nperbin) {
      tol.int <- 0.001
      int <- define.intervals(x = self$get.sVar(name.sVar), nbins = nbins, 
                              bin_bymass = bin_bymass, bin_bydhist = bin_bydhist, max_nperbin = max_nperbin)
      diffvec <- diff(int)
      if ( sum(abs(diffvec) < tol.int) > 0 ) {
        if (gvars$verbose) {
          message("No. of categories for " %+% name.sVar %+% " was collapsed from " %+%
                (length(int)-1) %+% " to " %+% (length(unique(int))-1) %+% " due to too few obs.")
          print("old intervals: "); print(int)
        }
        # int <- unique(int); print("new intervals: "); print(int)
        # Reason: Just taking unique interval values is insufficient
        # Instead need to drop all intervals that are "too close" to each other based on some tol value
        # remove all intervals (a,b) where |b-a| < tol.int, but always keep the very first interval (int[1])
        int <- c(int[1], int[2:length(int)][abs(diffvec) >= tol.int])
        if (gvars$verbose) print("new intervals: "); print(as.numeric(int))
      }
      return(int)
    },
    detect.cat.sVar.levels = function(name.sVar) {
      levels <- sort(unique(self$get.sVar(name.sVar)))
      return(levels)
    },
    # return the bin widths vector for the discretized continuous sVar (private$.ord.sVar):
    get.sVar.bw = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      intrvls.width <- diff(intervals)
      intrvls.width[intrvls.width <= gvars$tolerr] <- 1
      ord.sVar_bw <- intrvls.width[self$ord.sVar]
      return(ord.sVar_bw)
    },
   # return the bin widths differences vector for the discretized continuous sVar (self$ord.sVar):
    get.sVar.bwdiff = function(name.sVar, intervals) {
      if (!(self$active.bin.sVar %in% name.sVar)) stop("current discretized sVar name doesn't match name.sVar in get.sVar.bin.widths()")
      if (is.null(self$ord.sVar)) stop("sVar hasn't been discretized yet")
      ord.sVar_leftint <- intervals[self$ord.sVar]
      diff_bw <- self$get.sVar(name.sVar) - ord.sVar_leftint
      return(diff_bw)
    },
        
    # --------------------------------- binirize.sVar --------------------------------------
    # Purpose: Create a matrix of dummy bin indicators for categorical/continuous sVar
    # --------------------------------------------------------------------------------------
    binirize.sVar = function(name.sVar, ...) {
      private$.active.bin.sVar <- name.sVar
      if (self$is.sVar.cont(name.sVar)) {
        private$binirize.cont.sVar(name.sVar, ...)
      } else if (self$is.sVar.cat(name.sVar)) {
          if (self$is.sVar.CENS(name.sVar)) {
            private$binirize.cat.sVar(name.sVar, ...)
          } else {
            stop("...can only call $binirize.sVar for continuous or categorical sVars...")
          }
        }
    }, 
      
    # ---------------------------------- norm.cont.sVars -----------------------------------
    # Purpose: Normalize continuous sVars (** This could be memory-costly **)
    # --------------------------------------------------------------------------------------
    norm.cont.sVars = function() {
      names.c.sVar <- self$names.c.sVar
      if (length(names.c.sVar) == 0L) return(invisible(self))
      if (self$norm.c.sVars && (length(names.c.sVar) > 0)) {
        for (name.c.sVar in names.c.sVar) {
          self$dat.sVar[, name.c.sVar] <- normalize_sVar(self$dat.sVar[, name.c.sVar])
        }
      }
      invisible(self)
    },
      
    # ------------------------------------ fixmiss_sVar ------------------------------------
    # Purpose: Replace all missing (NA) values with a default integer (0)
    # --------------------------------------------------------------------------------------  
     fixmiss_sVar = function() {
      if (is.matrix(self$dat.sVar) | is.data.frame(self$dat.sVar)) { 
        private$fixmiss_sVar_mat()
      # } else if (is.data.table(self$dat.sVar)) { 
      #   private$fixmiss_sVar_DT()
      } else {
        stop("self$dat.sVar is of unrecognized class")
      }
    },
      
    # -------------------------------------- make.dat.sVar ---------------------------------------- 
    # Purpose: return a data.frames that combines all covariates (e.g., W, A, E);
    # --------------------------------------------------------------------------------------------- 
    # When !is.null(f.g_fun) create p new Odata's (n obs at a time) and replace A under g0 in Odata with A^* under g.star;
    # return a combined data.frame without saving individuals 
    # When is.null(f.g_fun), returns combined cbind(W, A) for observed Odata
    make.dat.sVar = function(p = 1, f.g_fun = NULL, regform = NULL)  {
      assert_that(is.count(p))
      self$p <- p
      nobs <- self$ndat.sVar
      WE.nodes <- self$nodes$WEnodes
      # set df.sVar to observed data with necessary nodes if g.fun is.null:
      if (is.null(f.g_fun)) {
        # assigning all covariates as data.frames:
        temp_dat.sVar <- self$dat.sVar  
        df.sVar <- temp_dat.sVar[, c(WE.nodes, self$nodes$Anodes)]
        # need to sample A under f.g_fun (gstar or known g0)
      } else {
        if (is.null(self$nodes$Anodes)) 
            stop("Anode(s) were not appropriately specified and is null; can't replace observed Anode with that sampled under g_star")
        df.sVar <- as.data.frame(matrix(nrow = (nobs * p), ncol = length(c(WE.nodes, self$nodes$Anodes))))  # pre-allocate result matx sVar
        temp_dat.sVar <- self$dat.sVar  
        WE_dat <- temp_dat.sVar[, WE.nodes]
        for (i in seq_len(p)) { 
          # if Anode is continuous, just call f.gen.probA.star:
          A.gstar <- f.gen.A.star(data = self$dat.sVar, f.g_fun = f.g_fun)
          # Replace A under g0 in Odata with A^* under g.star： Assigning the (W, E, A.gstar) to one output data matrix 
          df.sVar[((i - 1) * nobs + 1):(nobs * i), ] <-  cbind(WE_dat, A.gstar)[, ]
        }
        self$obs.wts <- rep(self$obs.wts, p)
      }
      colnames(df.sVar) <- c(WE.nodes, self$nodes$Anodes)  
      if (!is.null(regform) && length(setdiff(attributes(terms(regform))$term.labels, colnames(df.sVar))) > 0) {
        ExtraDat <- as.data.frame(model.matrix(as.formula(regform), data = df.sVar))
        df.sVar <- cbind(df.sVar, ExtraDat[, setdiff(names(ExtraDat), c(names(df.sVar), "(Intercept)")), drop = FALSE])
      }      
      self$dat.sVar <- as.data.frame(df.sVar)  # Caution: don't use data.frame(df.sVar) since it may change W2:A to W2.A
      invisible(self)
    }
  ),

  active = list(
    names.sVar = function() { colnames(self$dat.sVar) },
    names.c.sVar = function() { names(self$type.sVar[self$type.sVar %in% gvars$sVartypes$cont]) }, # names of cont sVars
    ncols.sVar = function() { length(self$names.sVar) },
    nobs = function() { nrow(self$dat.sVar) },
    dat.sVar = function(dat.sVar) {
      if (missing(dat.sVar)) {
        return(private$.mat.sVar)
      } else {
        assert_that(is.matrix(dat.sVar) | is.data.frame(dat.sVar)) # | is.data.table(dat.sVar) 
        private$.mat.sVar <- dat.sVar
      }
    },
    dat.bin.sVar = function(dat.bin.sVar) {
      if (missing(dat.bin.sVar)) {
        return(private$.mat.bin.sVar)
      } else {
        assert_that(is.matrix(dat.bin.sVar))
        private$.mat.bin.sVar <- dat.bin.sVar
      }
    },
    active.bin.sVar = function() { private$.active.bin.sVar },
    
    emptydat.sVar = function() { private$.mat.sVar <- NULL },  # wipe out mat.sVar
    emptydat.bin.sVar = function() {  # wipe out binirized mat.sVar:
      private$.mat.bin.sVar <- NULL
      private$.active.bin.sVar <- NULL
    },
    noNA.Ynodevals = function(noNA.Yvals) {
      if (missing(noNA.Yvals)) {
        return(private$.protected.YnodeVals)
      } else {
        private$.protected.YnodeVals <- noNA.Yvals
      }
    },
    nodes = function(nodes) {
      if (missing(nodes)) {
        return(private$.nodes)
      } else {
        assert_that(is.list(nodes))
        private$.nodes <- nodes
      }
    },    
    
    type.sVar = function() { private$.type.sVar }
  ),
          
  private = list(
    .nodes = list(),              # names of the nodes in the data (Anode, Wnode, etc..)
    .protected.YnodeVals = NULL,  # Actual observed values of the arbitrary outcome (Ynode in [0, 1]), along with deterministic vals
    .mat.sVar = NULL,             # pointer to data frame/matrix object storing the entire dataset (including all summaries sVars)
    .active.bin.sVar = NULL,      # Name of active binarized cont sVar, changes as fit/predict is called (bin indicators are temp. stored in mat.bin.sVar)
    .mat.bin.sVar = NULL,         # Temporary storage mat for bin indicators on currently binarized continous sVar (from private$.active.bin.sVar)
    .type.sVar = NULL,            # Named list with sVar types: list(names.sVar[i] = "binary"/"categor"/"contin"), can be overridden
    
    # Replace all missing (NA) values with a default integer (0) for matrix
    fixmiss_sVar_mat = function() {
      self$dat.sVar[gvars$misfun(self$dat.sVar)] <- gvars$misXreplace
      invisible(self)
    },
    # Replace all missing (NA) values with a default integer (0) for data frame/matrix/data.table
    fixmiss_sVar_DT = function() {
      # see http://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
      dat.sVar <- self$dat.sVar
      for (j in names(dat.sVar))
        set(dat.sVar, which(gvars$misfun(dat.sVar[[j]])), j , gvars$misXreplace)
      invisible(self)
    },
    # create a vector of ordinal (categorical) vars out of cont. sVar vector:
    discretize.sVar = function(name.sVar, intervals) {
      self$ord.sVar <- make.ordinal(x = self$get.sVar(name.sVar), intervals = intervals)
      invisible(self$ord.sVar)
    },
    # Create a matrix of bin indicators for continuous sVar:
    binirize.cont.sVar = function(name.sVar, intervals, nbins, bin.nms) {
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = private$discretize.sVar(name.sVar, intervals), nbins = nbins, bin.nms = bin.nms)
      invisible(self$dat.bin.sVar)
    },
    # Create a matrix of bin indicators for categorical sVar:
    binirize.cat.sVar = function(name.sVar, levels) {
      nbins <- length(levels)
      bin.nms <- self$bin.nms.sVar(name.sVar, nbins)
      self$dat.bin.sVar <- make.bins_mtx_1(x.ordinal = self$get.sVar(name.sVar), nbins = nbins, bin.nms = bin.nms, levels = levels)
      invisible(self$dat.bin.sVar)
    }
  )
)
