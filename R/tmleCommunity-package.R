#' Targeted Maximum Likelihood Estimation for Community-level Data
#' Finishing...
#' Finishing...
#' Finishing...

#' @section Documentation:
#' \itemize{
#' \item To see the package vignette use: \code{vignette("tmleCommunity_vignette", package="tmleCommunity")}
#' \item To see all available package documentation use: \code{help(package = 'tmleCommunity')}
#' }
#'
#' @section Routines:
#' The following routines will be generally invoked, in the same order as presented below.
#' \describe{
#'
#' Finishing...
#' }
#'
#' @section Datasets:
#' To learn more about the type of data input required by \code{\link{tmlenet}}, see the following example datasets:
#' \itemize{
#'   \item \code{\link{sampleDat_iidBinABinY}}
#'   \item \code{\link{sampleDat_iidBinAContY}}
#'   \item \code{\link{sampleDat_iidcontABinY}}
#'   \item \code{\link{sampleDat_iidcontAContY}}
#' }
#'
#' @section Updates:
#' Check for updates and report bugs at \url{https://github.com/chizhangucb/tmleCommunity}.
#'
#' @docType package
#' @name tmleCommunity-package
#'
NULL

#' An example of a continuous exposure with a continuous outcome.
#'
#' Simulated dataset containing measured i.i.d. baseline covariates (\code{W1}, \code{W2}, \code{W3} and \code{W4}), continuous  
#'  exposure (\code{A}) and continous outcome (\code{Y}). The 10,000 baseline covariates \code{W1}, \code{W2}, \code{W3} and \code{W4}
#'  were sampled as i.i.d., while the exposure value of \code{A} for each observation \code{i} was sampled conditionally on the value
#'  of \code{i}'s four baseline covariates, Similarly, the continuous outcome \code{Y} for each observation was generated conditionally 
#'  on \code{i}'s exposure and baseline covariates values in (\code{W1[i]},\code{W2[i]}, \code{W3[i]}, \code{W4[i]}, \code{A[i]}).
#'  Individual variables are described below.
#'
#' @format A data frame with 10,000 independent observations (rows) and 6 variables:
#' \describe{
#'   \item{W1}{binary baseline covariate with \code{P(W1 = 1) = 0.5}}
#'   \item{W2}{binary baseline covariate with \code{P(W1 = 1) = 0.3}}
#'   \item{W3}{continuous normal baseline covariate with mean = 0 and \eqn{\mu} = 0.25}}
#'   \item{W4}{continuous uniform baseline covariate with min = 0 and max = 1}}
#'   \item{A}{continuous normal exposure that depends on unit's baseline covariate values in \code{W1}, \code{W2}, \code{W3}, \code{W4}}
#'   \item{Y}{continuous normal  outcome that depends on unit's baseline covariate values and exposure in \code{W1}, \code{W2}, 
#'   \code{W3}, \code{W4}, \code{A}.

#' }
#' @docType data
#' @keywords datasets
#' @name sampleDat_iidcontAContY
#' @usage data(sampleDat_iidcontAContY)
NULL
