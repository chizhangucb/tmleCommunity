#' Targeted Maximum Likelihood Estimation for Community-level Data
#'
#' Targeted Maximum Likelihood Estimation (TMLE) for the sample average community-based treatment specific means effects (and can be  
#'  extended to additive treatment effect) in community-independent data structures. In other words, it estimates the marginal treatment  
#'  effect of single-time point arbitrary interventions on a continuous or binary outcome in community-independent data, adjusting for
#'  both community-level and individual-level baseline covariates. The package also provides Inverse-Probability-of-Treatment-Weighted 
#'  estimator (IPTW) and parametric G-computation formula estimator (GCOMP). The statistical inference (Standard errors, t statistc, 
#'  p-value and confidence intervals) of both TMLE and IPTW are based on the corresponding influence curve, respectively. Optional 
#‘  data-adaptive estimation of exposure and outcome mechanisms using the SuperLearner package and h2o package (latter for a large 
#'  dataset) is strongly recommended, especially when the outcome mechanism and treatment mechnism are unknown. Besides, it allows 
#'  for panel data transformation, such as with random effects and fixed effects. 
#' 
#' The input dataset should be made up of rows of community-specific and individual-specific observations, for community \eqn{j}, each  
#'  row \eqn{i} includes random variables \eqn{(W_{i,j}, E_{j}, A_{j}, Y_{i,j})}, where \eqn{E_j} represents a vector of community 
#'  \eqn{j}'s community-level (environmental) baseline covariates (individuals within the same community share the same values of 
#'  \eqn{E_j}), \eqn{W_{i,j}} represents a vector of individual \eqn{i}'s individual-level baseline covariates, \eqn{A_j} is the 
#'  exposure(s) (can be univariate or multivariate, can be binary, categorical or continuous) assigned or naturally occurred in  
#'  community \eqn{j} (individuals within the same community receive the same value of \eqn{A_j}) and \eqn{Y_{i,j}} is \eqn{i}'s 
#'  outcome (either binary or continuous). Each individual's baseline covariates \eqn{(W_{i,j}} depends on the environmental 
#'  baseline covariates \eqn{E_j} of the community \eqn{j} to which \eqn{i} belongs to. Similarly, each community's exposure 
#'  \eqn{A_j} depends on its community-level baseline covariates \eqn{E_j} and individual-level baseline covariates of all 
#'  individuals belonging to community \eqn{j} (all \eqn{W_{i,j}} such that \eqn{i} belongs to \eqn{j}). Besides, each outcome 
#'  \eqn{Y_{i,j}} could be affected by its baseline community and individual-level covariates \eqn{(E_j, W_{i,j})} and the baseline
#'  covariates of other individuals within the same community \eqn{(W_{s,j}: s\neq i, s\in j)}, together with its community-based
#'  intervention \eqn{A_j}. We note that the input data with no hierarchical structure (i.e., no communities and only individuals)
#'  is a special case of the hierarchical data since it simply treats \eqn{E_j} as \code{NULL}. 
#' 
#'  There are currently three approaches that can be used in hierarchical data analysis. The first community-level TMLE is developed 
#'  under a non-parametric causal model that allows for arbitrary interactions between individuals within a community. It estimates  
#'  the community-level causal effect by aggregating data at a community-level and treating community rather than the individual as 
#'  the unit of analysis (i.e., both community-level outcome and treatment mechanisms). The second individual-level TMLE is developed 
#'  under the submodel of the causal model in the first approach, incoporating knowledge of the dependence structure between 
#'  individual within communities (i.e., both individual-level outcome and treatmnet mechanisms). The third stratified TMLE fits a 
#'  separate outcome (exposure) mechanism for each community, and then combine those estimates into a (user-specific) average 
#'  (Default to be community size-weighed). Note that the stratified TMLE naturally controls for the community-level observed  
#'  covariates and unobserved factors. Namely, there is no \eqn{E} in the regressors for both outcome and treatment mechanisms.  
#' 
# @section Documentation:
# \itemize{
# \item To see the package vignette use: \code{vignette("tmleCommunity_vignette", package="tmleCommunity")}
# \item To see all available package documentation use: \code{help(package = 'tmleCommunity')}
# }
#'
#' @section References:
#' \enumerate{
#'  \item Balzer L. B., Zheng W., van der Laan M. J., Petersen M. L. and the SEARCH Collaboration (2017). A New Approach to 
#'    Hierarchical Data Analysis: Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference.
#'    ArXiv e-prints. 1706.02675.
#'  \item Mu\eqn{\~n}oz, I. D. and van der Laan, M. (2012). Population Intervention Causal Effects Based on Stochastic Interventions.
#'    Biometrics, 68(2):541-549.
#'  \item Sofrygin, O. and van der Laan, M. J. (2015). tmlenet: Targeted Maximum Likelihood Estimation for Network Data. 
#'    R package version 0.1.0.
#'  \item van der Laan, M. (2014). Causal Inference for a Population of Causally Connected Units. Journal of Causal Inference, 2(1)
#'  \item van der Laan, Mark J. and Gruber, Susan (2011). "Targeted Minimum Loss Based Estimation of an Intervention Specific 
#'    Mean Outcome". U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290. 
#'    http://biostats.bepress.com/ucbbiostat/paper290
#'  \item van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference for Observational and 
#'    Experimental Data" New York: Springer, 2011.
#'  \item Yves Croissant, Giovanni Millo (2008). Panel Data Econometrics in R: The plm Package. Journal of Statistical
#'.   Software 27(2). URL http://www.jstatsoft.org/v27/i02/.
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
#' 
#' To learn more about the type of data input required by \code{\link{tmleCommunity}}, see the following example datasets:
#' \itemize{
#'   \item \code{\link{comSample.wmT.bA.bY_list}}
#'   \item \code{\link{indSample.iid.cA.cY_list}}
#'   \item \code{\link{indSample.ind.bA.bY.rareJ1_list}}
#' }
#' There are a few other simulated datasets that can be used to test functions in \pkg{tmleCommunity}:
#' \itemize{
#'   \item \code{comSample.wmT.cA.cY_list}
#'   \item \code{indSample.iid.cA.bY_list}
#'   \item \code{\link{indSample.ind.bA.bY.rareJ2_list}}
#' }
#' @section Updates:
#' Check for updates and report bugs at \url{https://github.com/chizhangucb/tmleCommunity}.
#'
#' @docType package
#' @name tmleCommunity-package
#'
NULL

#' An Example of a Hierarchical Data Containing a Cluster-Based Binary Exposure with a Individual-Level Binary Outcome.
#'
#' Simulated hierarchical dataset containing 1000 independent communities, each (community \eqn{j}) containing \eqn{n_j} (non-fixed) 
#'  number of individuals where \eqn{n_j} is drawn from a normal with mean 50 and standard deviation 10 and round to the nearest 
#'  integer. Each row (observation) includes 2 measured community-level baseline covariates (\code{E1, E2}), 3 dependent   
#'  individual-level baseline covariates (\code{W1, W2, W3}), 1 dependent bianry exposure (\code{A}) and 1 dependent binary outcoem 
#'  (\code{Y}), along with one unique community identifier (\code{id}). The community-level baseline covariates (\code{E1, E2}) 
#'  were sampled as i.i.d across all communities, while the individual-level baseline covariates (\code{W1, W2, W3}) for each 
#'  individual \eqn{i} within communty \eqn{j} was generated conditionally on the values of \eqn{j}'s community-level baseline 
#'  covariates (\code{E1[j], E2[j]}). Then the community-level exposure (\code{A}) for each community \eqn{j} was sampled 
#'  conditionally on the value of \eqn{j}'s community-level baseline covariates (\code{E1[j], E2[j]}), together with all 
#'  invididuals' baseline covariates (\code{W1[i], W2[i], W3[i]}) within community \eqn{j} where \eqn{i=1,..,n_j}. Similary, 
#'  the individual-level binary outcome \code{Y} for each individual \eqn{i} within communty \eqn{j} was sampled conditionally 
#'  covariates and exposure (\code{E1[j], E2[j], A[j]}), as well as the value of individual \eqn{i}'s baseline covariates 
#'  on the value of community \eqn{j}'s baseline (\code{W1[i]}, \code{W2[i]}, \code{W3[i]}). The following section provides more 
#'  details regarding individual variables in simulated data.  R code for simulation of this dataset is at
#'  \url{https://github.com/chizhangucb/tmleCommunity/blob/master/tests/dataGeneration/get.cluster.dat.Abin.R}
#'
#' @format A data frame with 1000 independent communities, each containing around 50 individuals (in total 50,457 observations 
#'  (rows)), and 8 variables (columns):
#' \describe{
#'   \item{id}{integer (unique) community identifier from 1 to 1000, identical within the same community}
#'   \item{E1}{continuous uniform community-level baseline covariate with \code{min=0} and \code{max=1} (independent and identical
#'     across all individuals in the same community)}
#'   \item{E2}{discrete uniform community-level baseline covariate with 5 elements (0, 0.2, 0.4, 0.8, 1) (independent and identical
#'     across all individuals in the same community)}
#'   \item{W1}{binary individual-level baseline covariate that depends on the values of community-level baseline covaries (\code{E1,E2})}
#'   \item{W2}{continuous individual-level baseline covariate, together with \code{W3}, are drawn from a bivariate normal distribution
#'     with correlation 0.6, depending on the values of community's baseline covaries (\code{E1, E2})}
#'   \item{W3}{continuous normal individual-level baseline covariate, correlated with \code{W2}, see details in above}
#'   \item{A}{binary exposure that depends on community's baseline covariate values in \code{(E1, E2)}, and the mean of all individuals'
#'     baseline covariates \code{W1} within the same community}
#'   \item{Y}{binary outcome that depends on community's baseline covariate and exposure values in (\code{E1}, \code{E2}, \code{A}), 
#'     and all individuals' baseline covariate values in \code{(W2, W3)}}
#' }
#' @docType data
#' @keywords datasets
#' @name comSample.wmT.bA.bY_list
#' @usage data(comSample.wmT.bA.bY_list)
#'
#' @examples
#' data(comSample.wmT.bA.bY_list)
#' comSample.wmT.bA.bY <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
#' head(comSample.wmT.bA.bY)
#' comSample.wmT.bA.bY_list$psi0.Y  # 0.103716, True ATE
#' # summarize the number of individuals within each community
#' head(table(comSample.wmT.bA.bY$id))  
NULL

#' An Example of a Non-Hierarchical Data Containing a Continuous Exposure with a Continuous Outcome.
#'
#' Simulated (non-hierarchical) dataset containing 10,000 i.i.d. observations, with each row \code{i} consisting of measured baseline 
#'  covariates (\code{W1}, \code{W2}, \code{W3} and \code{W4}), continuous exposure (\code{A}) and continous outcome (\code{Y}). 
#'  The baseline covariates \code{W1}, \code{W2}, \code{W3} and \code{W4} were sampled as i.i.d., while the value of exposure \code{A} 
#'  for each observation \code{i} was drawn conditionally on the value of \code{i}'s four baseline covariates. Besides, the continuous
#'  outcome \code{Y} for each observation depends on \code{i}'s baseline covariates and exposure values in (\code{W1[i]},\code{W2[i]},
#'  \code{W3[i]}, \code{W4[i]}, \code{A[i]}). The following section provides more details regarding individual variables in simulated 
#'  data. R code for simulation of this dataset is at 
#'  \url{https://github.com/chizhangucb/tmleCommunity/blob/master/tests/dataGeneration/get.iid.dat.Acont.R}
#'
#' @format A data frame with 10,000 independent observations (rows) and 6 variables:
#' \describe{
#'   \item{W1}{binary baseline covariate with \eqn{P(W1=1) = 0.5}}
#'   \item{W2}{binary baseline covariate with \eqn{P(W2=1) = 0.3}}
#'   \item{W3}{continuous normal baseline covariate with \eqn{\mu} = 0 and \eqn{\sigma} = 0.25}
#'   \item{W4}{continuous uniform baseline covariate with \code{min=0} and \code{max=1}}
#'   \item{A}{continuous normal exposure where its mean depends on individual's baseline covariate values in \code{(W1, W2, W3, W4)}}
#      \code{W2}, \code{W3}, \code{W4}}
#'   \item{Y}{continuous normal outcome where its mean depends on individual's baseline covariate and exposure values in (\code{W1}, 
#'     \code{W2}, \code{W3}, \code{W4}, \code{A})}
#' }
#' @docType data
#' @keywords datasets
#' @name indSample.iid.cA.cY_list
#' @usage data(indSample.iid.cA.cY_list)
#' 
#' @examples
#' data(indSample.iid.cA.cY_list)
#' indSample.iid.cA.cY <- indSample.iid.cA.cY_list$indSample.iid.cA.cY
#' # True mean of outcome under intervention g0
#' psi0.Y <- indSample.iid.cA.cY_list$psi0.Y  
#' # True mean of outcoem under stochastic intervention gstar
#' psi0.Ygstar <- indSample.iid.cA.cY_list$psi0.Ygstar  
#' # truncated bound used in sampling A* under gstar (in data generating mechanism)
#' sampleDat_iidcontAContY$truncBD  
#' # shift value used in sampling A* under gstar 
#' sampleDat_iidcontAContY$shift.val
NULL

#' An Example of a Non-Hierarchical Data Containing a Binary Exposure with a Rare Binary Outcome (Independent Case-Control)
#'
#' Simulated (non-hierarchical) dataset containing 2,000 i.i.d. observations, with each row \code{i} consisting of 4 measured baseline 
#'  covariates (\code{W1}, \code{W2}, \code{W3} and \code{W4}), 1 binary exposure (\code{A}) and 1 binary outcome (\code{Y}) that
#'  defines case or control status. The baseline covariates \code{W1}, \code{W2}, \code{W3} and \code{W4} were sampled as i.i.d., 
#'  while the exposure \code{A} for each observation \code{i} depends on \code{i}'s four baseline covariates. Similarly, the outcome
#'  \code{Y} for each observation depends on \code{i}'s baseline covariates and exposure values. Moreover, we can also describe the 
#'  case-control design as first sampling \eqn{1} case \eqn{(W_1^1, W_2^1, W_3^1, W_4^1, A^1)} from the conditional distribution of 
#'  \eqn{(W_1, W_2, W_3, W_4, A)}, given Y = 1. One then samples \eqn{J} controls \eqn{(W_1^{0,j}, W_2^{0,j}, W_3^{0,j}, 
#'  W_4^{0,j}, A^{0,j})} from \eqn{(W_1, W_2, W_3, W_4, A)}, given Y = 0, \eqn{j=1,...,J}. Thus, the cluster containing one case 
#'  and \code{J} controls is considered the experimental unit. Finally one gets \eqn{nC} cases and \eqn{nCo} controls with 
#'  \eqn{J=nC/nCo}, where \eqn{J} can be used effectively in observation weights. The following section provides more details
#'  regarding individual variables in simulated data.
#'
#' @format A data frame with 2,000 independent observations (rows), containing 1000 cases and 1000 controls, and 6 variables, :
#' \describe{
#'   \item{W1}{continuous uniform baseline covariate with \code{min=0} and \code{max=1}}
#'   \item{W2}{continuous normal baseline covariate with \eqn{\mu} = 0 and \eqn{\sigma} = 0.3}
#'   \item{W3}{binary baseline covariate with \eqn{P(W2=1) = 0.5}}
#'   \item{W4}{binary baseline covariate with \eqn{P(W2=1) = 0.5}}
#'   \item{A}{binary exposure that depends on baseline covariate values in \code{(W1, W2, W3, W4)}}
#'   \item{Y}{binary outcome that depends on baseline covariate and exposure values in (\code{W1, W2, W3, W4, A})}
#'     
#' }
#' @docType data
#' @keywords datasets
#' @name indSample.ind.bA.bY.rareJ1_list
#' @usage data(indSample.ind.bA.bY.rareJ1_list)
#'
#' @examples
#' data(indSample.ind.bA.bY.rareJ1_list)
#' indSample.ind.bA.bY.rareJ1 <- indSample.ind.bA.bY.rareJ1_list$indSample.ind.bA.bY.rareJ1
#' head(indSample.ind.bA.bY.rareJ1_list$obs.wt.J1)  # Assigned weights to each observations
#' indSample.ind.bA.bY.rareJ1_list$q0  # 0.013579 True prevalence probability
#' indSample.ind.bA.bY.rareJ1_list$psi0.Y  # 0.012662 True ATE
#' indSample.ind.bA.bY.rareJ1_list$J  # 1 The ratio of number of controls to cases
NULL
