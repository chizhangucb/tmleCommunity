# tmleCommunity 
Targeted Maximum Likelihood Estimation for Hierarchical Data

## Summary

The `tmleCommunity` package performs targeted minimum loss-based estimation (TMLE) of the average causal effect of community-based intervention(s) at a single time point on an individual-based outcome of interest. It provides three approaches to analyze hierarchical data: community-level TMLE, inidividual-level TMLE and stratified TMLE. Implementations of the inverse-probability-of-treatment-weighting (IPTW) and the G-computation formula (GCOMP) are also available for each approach. The user-supplied arbitrary intervention can be either binary, categorical or continuous, also supporting univariate and multivariate setting. 

## Details

As a double-robust and asymptotically efficient substitution estimator that respects global constraints of the statistical model, targeted maximum likelihood (or minimum loss-based) estimation (TMLE) provides asymptotically valid statistical inference, with potential reduction in bias and gain in efficiency. The development of the `tmleCommunity` package for R was motivated by the increasing demand of a user-friendly tool to estimate the impact of community-based arbitrary exposures in community-independent data structures with a semi-parametric efficient estimator. Besides, the esimation results of TMLE, IPTW and GCOMP, the statistical inference (Standard errors, t statistc, p-value and confidence intervals) of both TMLE and IPTW are provided based on the corresponding influence curve, respectively. Optional data-adaptive estimation of exposure and outcome mechanisms using the `SuperLearner` package and `h2o` package (latter for a large dataset) is strongly recommended,

`tmleCommunity` is under active development so please submit any bug reports or feature requests to the [issue queue](https://github.com/chizhangucb/tmleCommunity/issues), or email Chi directly.

## Installation

### Github
To install the development version of tmleCommunity (requires the devtools package):

```{R install, eval=F}
# Install devtools if necessary:
# install.packages("devtools")
devtools::install_github("chizhangucb/tmleCommunity")
```

Alternatively, you could download the entire packge to the local path (e.g., Desktop) by either clicking the (green) download button on the this page, or cloning the repo through terminal via the following code

```console
git clone https://github.com/chizhangucb/tmleCommunity
```

Then open RStudio and run the following code 

```{R Load, eval=F}
# set the working directory to the directory where tmleCommunity pacakge is stored
setwd("some_path/tmleCommunity")

# 1. If you only want to use the package instead of installing it in R library, use 
devtools::load_all()

# 2. If you want to install it, then uze
devtools::install()
library(tmleCommunity)
```

### CRAN

Forthcoming Summer 2018

## Documentation

Once the package is installed, please refer to the help file `?'tmleCommunity-package'` and `tmleCommunity` function documentation for details and examples.

```{R Documentation, eval=F}
?'tmleCommunity-package'
?tmleCommunity
```

## Example

We will use the sample dataset (`E`=(`E1`,`E2`),`W`=(`W1`,`W2`,`W3`),`A`,`Y`) that come along with the package:

```{R Data, eval=F}
data(comSample.wmT.bA.bY_list)  # load the sample data 
comSample.wmT.bA.bY <- comSample.wmT.bA.bY_list$comSample.wmT.bA.bY
N <- NROW(comSample.wmT.bA.bY)
```

Estimating the additive treatment effect (ATE) for two deterministic interventions (`f_gstar1 = 1` vs `f_gstar2 = 0`) via community-level / individual-level analysis.

```{R GLM_analysis, eval=F}
# speed.glm using correctly specified Qform, hform.g0 and hform.gstar;
Qform.corr <- "Y ~ E1 + E2 + W2 + W3 + A" # correct Q form
gform.corr <- "A ~ E1 + E2 + W1"  # correct g

# Setting global options that may be used in tmleCommunity(), e.g., using speed.glm
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", maxNperBin = N)

# Community-level analysis without a pooled individual-level regression on outcome
tmleCom_wmT.bA.bY.1a_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = FALSE, 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)


# Community-level analysis with a pooled individual-level regression on outcome
tmleCom_wmT.bA.bY.1b_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = TRUE, 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.1b_sglm$ATE$estimates

# Individual-level analysis with both individual-level outcome and treatment mechanisms
tmleCom_wmT.bA.bY.2_sglm <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "individual_level", communityID = "id", 
                Qform = Qform.corr, hform.g0 = gform.corr, hform.gstar = gform.corr)
tmleCom_wmT.bA.bY.2_sglm$ATE$estimates
```

If you are uncertain about the model specification of exposure and outcome mechanisms, data-adaptive estimation methods may be a better choice than parametric models. 

```{R Data_adaptive_analysis, eval=F}
# SuperLearner for both outcome and treatment (clever covariate) regressions
# using all parent nodes (of Y and A) as regressors (respectively)
require("SuperLearner")
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", 
                maxNperBin = N, SL.library = c("SL.glm", "SL.step", "SL.bayesglm"))
tmleCom_wmT.bA.bY.2_SL <- 
  tmleCommunity(data = comSample.wmT.bA.bY, Ynode = "Y", Anodes = "A", 
                WEnodes = c("E1", "E2", "W1", "W2", "W3"), f_gstar1 = 1L, f_gstar2 = 0L,
                community.step = "community_level", communityID = "id", pooled.Q = TRUE, 
                Qform = NULL, hform.g0 = NULL, hform.gstar = NULL)
tmleCom_wmT.bA.bY.2_SL$ATE$estimates
```

## Authors

Chi Zhang, Oleg Sofrygin, Jennifer Ahern, M. J. van der Laan

## References

Balzer L. B., Zheng W., van der Laan M. J., Petersen M. L. and the SEARCH Collaboration (2017). A New Approach to Hierarchical Data Analysis: Targeted Maximum Likelihood Estimation of Cluster-Based Effects Under Interference. ArXiv e-prints. 1706.02675.

Muoz, I. D. and van der Laan, M. (2012). Population Intervention Causal Effects Based on Stochastic Interventions. Biometrics, 68(2):541-549.

Sofrygin, O. and van der Laan, M. J. (2015). tmlenet: Targeted Maximum Likelihood Estima- tion for Network Data. R package version 0.1.9. https://github.com/osofr/tmlenet

van der Laan, M. (2014). Causal Inference for a Population of Causally Connected Units. Journal of Causal Inference, 2(1)

van der Laan, Mark J. and Gruber, Susan (2011). "Targeted Minimum Loss Based Estimation of an Intervention Specific Mean Outcome". U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 290. http://biostats.bepress.com/ucbbiostat/paper290

van der Laan, Mark J. and Rose, Sherri, "Targeted Learning: Causal Inference for Observa- tional and Experimental Data" New York: Springer, 2011.


## News 

`tmleCommunity` 0.1.0

=====================

2018-08-14

`h2oEnsemble` has been removed from the list of required libraries since it's not in mainstream repositories.
