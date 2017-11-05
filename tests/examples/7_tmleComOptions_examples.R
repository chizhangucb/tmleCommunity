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

# library including "SL.glm", "SL.glmnet", "SL.gam", and "SL.randomForest"
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.gam", "SL.randomForest"), CVfolds = 10)

# Create custom learners (e.g., a sequence of learners with hyperparameter combinations)


tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = 10000,
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), CVfolds = 10)
                


tmleCom_Options(Qestimator = "SuperLearner", gestimator = "h2o__ensemble", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
                h2ometalearner = "h2o.deeplearning.wrapper", 
                h2olearner = c("h2o.gbm.wrapper", "h2o.randomForest.wrapper"))
