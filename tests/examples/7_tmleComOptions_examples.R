#***************************************************************************************
# Example 1: using different estimators in estimation of Q and g mechanisms
#***************************************************************************************
# 1.1 using speed.glm (and glm)
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm")
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "glm__glm")

# 1.2 using uperLearner
# library including "SL.glm", "SL.glmnet", "SL.ridge", and "SL.stepAIC"
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", CVfolds = 5,
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"))

# library including "SL.bayesglm", "SL.gam", and "SL.randomForest", and split to 10 CV folds
# require("gam"); require("randomForest")
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", CVfolds = 10,
                SL.library = c("SL.bayesglm", "SL.gam", "SL.randomForest"))

# Create glmnet wrappers with different alphas (the default value of alpha in SL.glmnet is 1)
create.SL.glmnet <- function(alpha = c(0.25, 0.50, 0.75)) {
  for(mm in seq(length(alpha))){
    eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', 
                            alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), 
         envir = .GlobalEnv)
  }
  invisible(TRUE)
}
create.SL.glmnet(seq(0, 1, length.out=3))  # 3 glmnet wrappers with alpha = 0, 0.5, 1
# Create custom randomForest learners (set ntree to 100 rather than the default of 500) 
create.SL.rf <- create.Learner("SL.randomForest", list(ntree = 100))
# Create a sequence of 3 customized KNN learners 
# set the number of nearest neighbors as 8 and 12 rather than the default of 10
create.SL.Knn <- create.Learner("SL.kernelKnn", detailed_names = T, tune = list(k = c(8, 12)))
SL.library <- c(grep("SL.glmnet.", as.vector(lsf.str()), value=TRUE), 
                create.SL.rf$names, create.SL.Knn$names)
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", 
                SL.library = SL.library, CVfolds = 5)            

# 1.3 using h2o.ensemble
# h2olearner including "h2o.glm.wrapper" and "h2o.randomForest.wrapper"
require("h2oEnsemble")
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", 
                CVfolds = 10, h2ometalearner = "h2o.glm.wrapper", 
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"))

# Create a sequence of customized h2o glm, randomForest and deeplearning wrappers 
h2o.glm.1 <- function(..., alpha = 1, prior = NULL) { 
  h2o.glm.wrapper(..., alpha = alpha, , prior=prior) 
}
h2o.glm.0.5 <- function(..., alpha = 0.5, prior = NULL) { 
  h2o.glm.wrapper(..., alpha = alpha, , prior=prior) 
}
h2o.randomForest.1 <- function(..., ntrees = 200, nbins = 50, seed = 1) {
  h2o.randomForest.wrapper(..., ntrees = ntrees, nbins = nbins, seed = seed)
}
h2o.deeplearning.1 <- function(..., hidden = c(500, 500), activation = "Rectifier", seed = 1) {
  h2o.deeplearning.wrapper(..., hidden = hidden, activation = activation, seed = seed)
}
h2olearner <- c("h2o.glm.1", "h2o.glm.0.5", "h2o.randomForest.1", 
                "h2o.deeplearning.1", "h2o.gbm.wrapper")
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble",
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
                h2ometalearner = "h2o.deeplearning.wrapper", h2olearner = h2olearner)

# using "h2o.deeplearning.wrapper" for h2ometalearner
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble",
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
                h2ometalearner = "h2o.deeplearning.wrapper", h2olearner = h2olearner)

#***************************************************************************************
# Example 2: Define the values of bin cutoffs for continuous outcome in different ways
# through three arguments - bin.method, nbins, maxNperBin 
#***************************************************************************************
# 2.1 using equal-length method
# discretize a continuous outcome variable into 10 bins, no more than 1000 obs in each bin 
tmleCom_Options(bin.method = "equal.len", nbins = 10, maxNperBin = 1000)

# 2.2 find a compromise between equal-mass and equal-length method
# discretize into 5 bins (default), and no more than 5000 obs in each bin
tmleCom_Options(bin.method = "dhist", nbins = 10, maxNperBin = 5000)

# 2.3 Default to use equal-mass method with 5 bins, no more than 500 obs in each bin
tmleCom_Options()
