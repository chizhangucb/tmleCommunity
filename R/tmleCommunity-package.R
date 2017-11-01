#' Targeted Maximum Likelihood Estimation for Community-level Data
#'
#' The \textbf{tmleCommunity} package implements the Targeted Maximum Likelihood Estimation (TMLE) for the sample average community-based  
#'  treatment specific means effects (and can be extended to ATE) in community-independent data structures. In other words, it estimates  
#'  the marginal treatment effect of single-time point arbitrary interventions on a continuous or binary outcome in community-independent 
#'  data, adjusting for both community-level and individual-level baseline covariates. The package also provides Inverse-Probability-of-
#'  Treatment-Weighted estimator (IPTW) and parametric G-computation formula estimator (GCOMP). The statistical inference (Standard errors, 
#'  p-value and confidence intervals) of both TMLE and IPTW are based on the corresponding influence curve, respectively. Optional data-adaptive 
#‘  estimation of exposure and outcome mechanisms using the SuperLearner package and h2o package (latter for a large dataset) is strongly 
#'  recommended, especially when the outcome mechanism and treatment mechnism are unknown. Besides, it allows for panel data transformation, 
#'  such as with random effects and fixed effects. 
#' 
#' The input dataset should be made up of rows of unit-specific observations, each row i includes variables (W_i, E_i, A_i, Y_i), where W_i 
#'  represents a vector of i’s individual-level baseline covariates, E_i represents a vector of i’s community-level baseline covariates 
#'  (observations within the same community usually have the same values of E_i), A_i is a vector of i’s interventions (can be univariate 
#'  or multivariate, can be binary, categorical or continuous), and Y_i is i’s outcome. For each community, individual exposure and outcome
#'  mechanisms will be estimated, then the ATE across all the communities is calculated as a user-specific average of all community-level 
#'  estimates (Default to size-weighted). Besides, each exposure A_i is a function of baseline covariates (W_i, E_i), and the outcome Y_i 
#'  is a function of both baseline and exposure covariates (W_i, E_i, A_i). 
#' 
#' 
#' @section Documentation:
#' \itemize{
#' \item To see the package vignette use: \code{vignette("tmleCommunity_vignette", package="tmleCommunity")}
#' \item To see all available package documentation use: \code{help(package = 'tmleCommunity')}
#' }
#'
#' @section Reference(s):
#' \enumerate{
#'   \item Muñoz, I. D. and van der Laan, M. (2012). Population intervention causal effects based on stochastic interventions. Biometrics, 
#'         68(2):541–549.
#'   \item Sofrygin, O. and van der Laan, M. J. (2015). tmlenet: Targeted Maximum Likelihood Estimation for Network Data. R package version 0.1.0.
#'   \item van der Laan, Mark J. and Gruber, Susan, "Targeted Minimum Loss Based Estimation of an Intervention Specific Mean Outcome" 
#'         (August 2011). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290. 
#'         http://biostats.bepress.com/ucbbiostat/paper290
#'   \item van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference for Observational and Experimental Data" New York: 
#'         Springer, 2011.
#' }
#'
# @section Routines:
# The following routines will be generally invoked, in the same order as presented below.
# \describe{
#
# Finishing...
# }
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
#'   \item{W3}{continuous normal baseline covariate with mean = 0 and \eqn{\mu} = 0.25}
#'   \item{W4}{continuous uniform baseline covariate with min = 0 and max = 1}
#'   \item{A}{continuous normal exposure that depends on unit's baseline covariate values in \code{W1}, \code{W2}, \code{W3}, \code{W4}}
#'   \item{Y}{continuous normal  outcome that depends on unit's baseline covariate values and exposure in \code{W1}, \code{W2}, 
#'   \code{W3}, \code{W4}, \code{A}}
#' }
#' @docType data
#' @keywords datasets
#' @name sampleDat_iidcontAContY
#' @usage data(sampleDat_iidcontAContY)
#' @examples
#'
#' data(sampleDat_iidcontAContY)
#' dat_iidcontAContY <- sampleDat_iidcontAContY$dat_iidcontAContY
#' psi0.Y <- sampleDat_iidcontAContY$psi0.Y
#' psi0.Ygstar <- sampleDat_iidcontAContY$psi0.Ygstar
#' sampleDat_iidcontAContY$truncBD
#' sampleDat_iidcontAContY$shift.val
NULL
