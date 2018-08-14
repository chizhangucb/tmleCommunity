library(testthat)
library(tmleCommunity)
library(devtools)
install_github("h2oai/h2o-3/h2o-r/ensemble/h2oEnsemble-package")
library(h2oEnsemble)

test_check("tmleCommunity")
