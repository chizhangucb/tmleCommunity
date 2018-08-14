# tmleCommunity 
Targeted Maximum Likelihood Estimation for Hierarchical Data

## Summary

The `tmleCommuniy` package performs targeted minimum loss-based estimation (TMLE) of the average causal effect of community-based intervention(s) at a single time point on an individual-based outcome of interest. It provides three approaches to analyze hierarchical data: community-level TMLE, inidividual-level TMLE and stratified TMLE. Implementations of the inverse-probability-of-treatment-weighting (IPTW) and the G-computation formula (GCOMP) are also available for each approach. The user-supplied arbitrary intervention can be either binary, categorical or continuous, also supporting univariate and multivariate setting. 

## Details

As a double-robust and asymptotically efficient substitution estimator that respects global constraints of the statistical model, targeted maximum likelihood (or minimum loss-based) estimation (TMLE) provides asymptotically valid statistical inference, with potential reduction in bias and gain in efficiency. The development of the `tmleCommunity` package for R was motivated by the increasing demand of a user-friendly tool to estimate the impact of community-based arbitrary exposures in community-independent data structures with a semi-parametric efficient estimator. Besides, the esimation results of TMLE, IPTW and GCOMP, the statistical inference (Standard errors, t statistc, p-value and confidence intervals) of both TMLE and IPTW are provided based on the corresponding influence curve, respectively. Optional data-adaptive estimation of exposure and outcome mechanisms using the `SuperLearner` package and `h2o` package (latter for a large dataset) is strongly recommended,

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
