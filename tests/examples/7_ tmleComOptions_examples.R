tmleCom_Options(Qestimator = "speedglm__glm", gestimator = "speedglm__glm", nbins = 20)

tmleCom_Options(Qestimator = "h2o__ensemble", gestimator = "h2o__ensemble", maxNperBin = 10000,
                h2olearner = c("h2o.glm.wrapper", "h2o.randomForest.wrapper"), CVfolds = 10)
                
tmleCom_Options(Qestimator = "SuperLearner", gestimator = "SuperLearner", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5)

tmleCom_Options(Qestimator = "SuperLearner", gestimator = "h2o__ensemble", maxNperBin = nrow(data),
                SL.library = c("SL.glm", "SL.glmnet", "SL.ridge", "SL.stepAIC"), CVfolds = 5,
                h2ometalearner = "h2o.deeplearning.wrapper", 
                h2olearner = c("h2o.gbm.wrapper", "h2o.randomForest.wrapper"))
