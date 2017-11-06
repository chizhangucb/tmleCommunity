#***************************************************************************************
# Example 1: using different estimators in estimation of Q and g mechanisms
#***************************************************************************************
# 1.1 speed.glm (and glm)
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm")
tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "glm__glm")

# 1.2 SuperLearner
# library including "SL.glm", "SL.glmnet", "SL.ridge", and "SL.stepAIC"
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5)

# library including "SL.bayesglm", "SL.gam", and "SL.randomForest", and split to 10 CV folds
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
                SL.library = c("SL.bayesglm", "SL.gam", "SL.randomForest"), CVfolds = 10)

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
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
                SL.library = SL.library, CVfolds = 5)                
                
tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = 10000,
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), CVfolds = 10)
                
# 1.3 SuperLearner

tmleCom_Options(Qestimator = "SuperLearner", gestimator = "h2o__ensemble", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
                h2ometalearner = "h2o.deeplearning.wrapper", 
                h2olearner = c("h2o.gbm.wrapper", "h2o.randomForest.wrapper"))
