# tmleCommunity - Targeted Maximum Likelihood Estimation for Hierarchical Data

The tmleCommuniy package performs targeted minimum loss-based estimation (TMLE) of the average causal effect of community-based intervention(s) at a single time point on an individual-based outcome of interest. It provides three approaches to analyze hierarchical data: community-level TMLE, inidividual-level TMLE and stratified TMLE. Implementations of the inverse-probability-of- treatment-weighting (IPTW) and the G-computation formula (GCOMP) are also available for each approach. The user-supplied arbitrary intervention can be either binary, categorical or continuous, also supporting univariate and multivariate setting. 

## Installation and Documentation

### Install from github
To install the development version of tmleCommunity (requires the devtools package):

devtools::install_github("chizhangucb/tmleCommunity")

### Install from local source
Alternatively, could download the entire packge to the local path (e.g., Desktop) by either clicking the (green) download button on the this page, or cloning the repo through terminal via the following code

git clone https://github.com/chizhangucb/tmleCommunity

Then open RStudio and set the working directory to the directory where tmleCommunity pacakge is stored, via 

setwd("some-path/tmleCommunity")

If you only want to use the package instead of installing it in R library, use 

devtools::load_all()

If you want to install it, then 

devtools::install()

library(tmleCommunity)

### Generate package document 
Once you set the working directory to the directory where tmleCommunity pacakge is stored, use 

devtools::document()
