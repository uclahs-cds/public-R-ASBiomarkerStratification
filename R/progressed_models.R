library(caret)
# library(randomForest)
library(rpart)
library(rpart.plot)
# library(C50)
library(ProstateCancer.ASBiomarkerSynergy)
library(pROC)
# Boosted Trees
library(gbm)

set.seed(13)

biodb <- default.load.data(onlyBiodb = TRUE);

# No training/test since using cross-validation
# inTrain <- createDataPartition(
#     y = biodb$ProgressedToTreatment,
#     p = .8,
#     list = FALSE
#     );
#
# str(inTrain)
#
# training <- biodb[inTrain, ]
# testing <- biodb[-inTrain, ]

biokey.variables <- c('Age',
                      'Race',
                      'Ethnicity',
                      'Weight',
                      'Height',
                      'BMI', 'MRIResult', 'MRILesions',
                      # 'BiopsyResult',
                      'ProstateVolume',
                      # 'Observation',
                      'PCA3',
                      'T2ERG', 'MiPSCancerRisk', 'MiPSHighGradeCancerRisk', 'PSAHyb',
                      'freePSA', 'p2PSA', 'PercentFreePSA', 'PHI', 'GeneticAncestry',
                      # 'PreviousGleason', 'StudyHighestGleason', 'RSIlesionGleason', # Only use the ISUP grade group over Gleason
                      # 'StudyHighestISUP', 'HighestPIRADS',
                      'PreviousISUP',
                      'GeneticRiskScore',
                      'TNFaAverage', 'TNFaSTD', # Tumor necrosis factor?
                      #'GeneticRiskCategory', Don't need since genetic risk score is continuous version of this
                      # 'GlobalScreeningArray', This is just an indicator if any of the follow are > 0
                      'GSAPositives', 'BRCAMutation',
                      'Mutation.BRCA1', 'Mutation.BRCA2', 'Mutation.ATM', 'Mutation.MLH1', 'Mutation.PMS2',
                      'RSInormalSignal',
                      'RSIlesionSignal',
                      'ADCnormalSignal', 'ADClesionSignal',
                      # 'RSIlesionPIRADS', 'RSIlesionCancer',  'RSIlesionUpgraded', 'RSIlesionISUP',
                      'PSADensity', 'PHIDensity');

gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                        n.trees = (1:30)*50,
                        shrinkage = c(0.001, 0.01, 0.1),
                        n.minobsinnode = 20)

custom.summary <- function (data, lev = NULL, model = NULL) {
    c(
        twoClassSummary(data, lev, model),
        prSummary(data, lev, model)
    )
}

ctrl.sampling <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = custom.summary,
    # https://topepo.github.io/caret/subsampling-for-class-imbalances.html
    sampling = "up"
)

y <- biodb$ProgressedToTreatment;
levels(y) <- c("no", "yes");

rpart.fit <- train(
    biodb[, biokey.variables],
    y,
    method = "rpart",
    metric = "F",
    trControl = ctrl.sampling,
    tuneLength = 30
)

rpart.plot(rpart.fit$finalModel);

c50.grid <- expand.grid(
    .winnow = c(TRUE, FALSE),
    .trials = c(1, 5, 10, 15, 20),
    .model = c('tree', 'rules')
)

c50.fit <- train(
    biodb[, biokey.variables],
    y,
    method = "C5.0",
    metric = "F",
    trControl = ctrl.sampling,
    tuneGrid = c50.grid
)

c50.params <- as.list(c50.fit$bestTune)
c50.params$x <- biodb[, biokey.variables]
c50.params$y <- biodb$ProgressedToTreatment

c50.bestfit <- do.call(C5.0, c50.params)
plot(c50.bestfit)

varImp(c50.bestfit)

# Remove RSI variables
# biokey.variables.no.RSI <- biokey.variables[!grepl("RSI", biokey.variables)]
#
# model.formula <- as.formula(paste('ProgressedToTreatment', paste0(biokey.variables, collapse = " + "), sep = "~"))
# model.formula.bx <- as.formula(paste('BiopsyUpgraded', paste0(biokey.variables, collapse = " + "), sep = "~"))
# model.formula.bx.no.RSI <- as.formula(paste('BiopsyUpgraded', paste0(biokey.variables.no.RSI, collapse = " + "), sep = "~"))
#
# mod <- rpart(model.formula, data = biodb)
# modC5 <- C5.0(model.formula, data = biodb)
#
# mod$variable.importance
# rpart.plot(mod)
# mod.preds <- predict(mod, type = "class")
# confusionMatrix(mod.preds, biodb$ProgressedToTreatment)
# confusionMatrix(predict(modC5, newdata = biodb, type = "class"), biodb$ProgressedToTreatment)
# summary(modC5)
#
# ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3, classProbs = TRUE)
#
# mod.bx.train <- train(biodb[, biokey.variables],
#                 biodb$BiopsyUpgraded,
#                 method = "rpart",
#                 trControl = ctrl)
# mod.bx <- rpart(model.formula.bx, data = biodb, cp = mod.bx.train$bestTune)
# rpart.plot(mod.bx)
# summary(mod.bx)
#
# mod.bx <- C5.0(model.formula.bx, data = biodb)
# mod.bx.no.RSI <- C5.0(model.formula.bx.no.RSI, data = biodb)
# summary(mod.bx)
# summary(mod.bx.no.RSI)
#
# mod.bx.preds <- predict(mod.bx, newdata = biodb, type = "class")
# confusionMatrix(mod.bx.preds[!is.na(biodb$BiopsyUpgraded)], biodb$BiopsyUpgraded[!is.na(biodb$BiopsyUpgraded)])
#
# mod.bx.no.RSI.preds <- predict(mod.bx.no.RSI, newdata = biodb, type = "class")
# confusionMatrix(mod.bx.no.RSI.preds[!is.na(biodb$BiopsyUpgraded)], biodb$BiopsyUpgraded[!is.na(biodb$BiopsyUpgraded)])
#
# ## Caret
# rpart.train <- train(biodb[, biokey.variables],
#                      y,
#                      method = "rpart",
#                      trControl = ctrl2,
#                      metric = "F",
#                      tuneLength = 30
#                      )
#
# mod.f <- rpart(model.formula, data = biodb, cp = rpart.train$bestTune)
# mod.f.preds <- predict(mod.f, type = "class")
# confusionMatrix(data = mod.f.preds, reference = biodb$ProgressedToTreatment)

gbm.fit.progress <- train(
    biodb[, biokey.variables],
    y,
    method = "gbm",
    metric = "F",
    trControl = ctrl.sampling,
    tuneGrid = gbmGrid,
    verbose = FALSE
)

summary(gbm.fit.progress)

### Compare rpart, C5.0, and GBM
resamps <- resamples(list(GBM = gbm.fit.progress,
                          C5.0 = c50.fit,
                          rpart = rpart.fit))

summary(resamps)

# Variable importance
# lapply(list(rpart.fit, c50.fit, gbm.fit.progress), varImp)
# do.call(cbind, lapply(list(rpart.fit, c50.fit, gbm.fit.progress), function(x) varImp(x)$importance))
