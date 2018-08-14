## Initial Submission to CRAN:

## R CMD check --as-cran results:
There were one ERROR:

* Package required but not available: 'h2oEnsemble'

This package provides important additional functionality to `tmleCommunity`. Source code for h2oEnsemble is available from \url{https://github.com/h2oai/h2o-3/tree/master/h2o-r/ensemble} 

There was one WARNING:

* checking CRAN incoming feasibility ... WARNING

Maintainer: 'Chi Zhang <chi.zhang@berkeley.edu>'

According to CRAN Maintainer Uwe Ligges, *This is just a note that reminds CRAN maintainers 
to check that the submission comes actually from his maintainer and not anybody else.*. 
Thus, it is safe to ignore such a message.

There was one NOTE:

* Possibly mis-spelled words in DESCRIPTION:
  GCOMP (17:63)
  IPTW (17:26)
  SuperLearner (22:41)
  TMLE (13:54, 15:90, 16:23, 16:43, 21:41)
  i's (26:28, 28:54)
  inidividual (16:5)
  j's (25:42)
  oEnsemble (22:60)
  
These are not mis-spelled. I manually checked all of these words.
They define specific estimators / terms that are known to the intended user.  
Caution: oEnsemble is not a word, it's supposed to be h2oEnsemble.




