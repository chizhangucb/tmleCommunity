#----------------------------------------------------------
# General utilities / Global Vars
#----------------------------------------------------------
is.integerish <- function (x) { is.integer(x) || (is.numeric(x) && all(x == as.integer(x))) }

iqr <- function(x) { return(diff(quantile(x,c(.25, .75),na.rm=T))) }  # interquartile range

## Original method code from https://github.com/cran/plm/blob/master/R/plm.R
## https://github.com/cran/plm/blob/master/R/pFormula.R

# Discarding missing outcomes (by default).

#------------------------------------
#' Transfer a panel dataset into a fixed-effect/ random-effect transformed data, using individual (and time) indexes
#'
#' \code{panelData.Trans} provides a wide variety of ways of data transformation for panel datasets. It allows users 
#'  to only apply transformation on regressors of interests, instead of on the entire dataset. See details in 
#'  \url{https://github.com/cran/plm/blob/master/R/plm.R}.
#' @param data A data frame (will be automatically transferred to panel data frame) or a panel data frame
#' @param yvar Column names in \code{data} of outcome variable (Currently only support univariate).    
#' @param xvar Column names in \code{data} of explanatory variables (Including \eqn{(A, W, E)}).
#' @param effect The effects introduced in the model, one of "individual", "time", "twoways" and "nested". Default to "individual".
#' @param model Model of estimation, one of "pooling" (pooled OLS), "within" (fixed effect), "between" (group mean), 
#'  "random"(random effect), "fd" (first differences) and "ht" (Hausman-Taylor estimator). Default to "within".
#' @param index A vector of two character strings which contains the names of the individual and of the time indices. 
#'  If only individual index is given, treat each observation within a unit as a time point.
#'  If no index is given, the first two columns will be automatically treated as individual and time indices, sequentially.
#' @param transY Logical. If \code{TRUE}, indicate the outcome variable \code{yvar} will also be tranformed. Default to \code{TRUE}.
#' @return 
#' \itemize{
#'   \item \code{newdata} Transformed panel data
#' }
#' @export
panelData.Trans <- function(data, yvar, xvar, effect = "individual", model = "within", index = NULL, transY = TRUE) {
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
  X <- plm:::model.matrix.pFormula(object = formula, data, model = model, effect = effect, rhs = 1, theta = NULL)
  if (ncol(X) == 0) stop("empty model")
  if (transY) {
    y <- plm:::pmodel.response(object = formula, data, model = model, effect = effect, theta = NULL)
  } else {
    y <- data[, yvar]
  }
  
  newdata <- as.data.frame(cbind(y, X))
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
  assert_that(is.integerish(getopt("max
