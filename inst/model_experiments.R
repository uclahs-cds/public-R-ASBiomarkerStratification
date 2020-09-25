library(caret)
library(ProstateCancer.ASBiomarkerSynergy);
library(rpart.plot);

biodb <- default.load.data(onlyBiodb = TRUE);

biokey.variables <- c('Age',
                      'Race',
                      'Hispanic',
                      'Weight',
                      'Height',
                      'BMI', 'MRIResult', 'MRILesions',
                      'BiopsyResult',
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
                      'GlobalScreeningArray', # This is just an indicator if any of the follow are > 0
                      'GSAPositives', 'BRCAMutation',
                      'Mutation.BRCA1', 'Mutation.BRCA2', 'Mutation.ATM', # 'Mutation.MLH1', 'Mutation.PMS2', Not enough data
                      'RSInormalSignal',
                      'RSIlesionSignal',
                      'ADCnormalSignal', 'ADClesionSignal',
                      # 'RSIlesionPIRADS', 'RSIlesionCancer',  'RSIlesionUpgraded', 'RSIlesionISUP',
                      'PSADensity', 'PHIDensity');

target <- "ProgressedToTreatment"
missing.target <- is.na(biodb[, target])
# Variable we want to predict
y.target <- biodb[!missing.target, target];

y <- y.target
levels(y) <- c('no', 'yes');

# Drop the missing target values
X <- biodb[!missing.target, biokey.variables]
X.impute.fit <- preProcess(X, method = "knnImpute", na.remove = TRUE)
X.imputed <- predict(X.impute.fit, X)

custom.summary <- function (data, lev = NULL, model = NULL) {
    c(
        defaultSummary(data, lev, model),
        twoClassSummary(data, lev, model),
        prSummary(data, lev, model)
    )
}

train.control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = custom.summary
)

train.control.up <- train.control
train.control.up$sampling <- 'up'

metric <-  'F'
rpart.fit <- train(
    X,
    y,
    method = 'rpart',
    metric = metric,
    trControl = train.control.up,
    tuneLength = 30
)

confusionMatrix(predict(rpart.fit), reference = y)

rforest.fit <- train(
    X.imputed[,apply(X.imputed, 2, function(x) !any(is.na(x)))],
    y,
    method = 'rf',
    metric = metric,
    trControl = train.control,
    tuneLength = 30
)


tune.grid <- expand.grid(
    nrounds = seq(from = 200, to = 1000, by = 50),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 3, 4, 5, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
)

prog.form <- as.formula(paste0("~ ", paste0(biokey.variables, collapse = " + ")))
X.dummy <- predict(dummyVars(prog.form, data = biodb, fullRank = TRUE), newdata = biodb)

# Convert ordered factors to numeric rather than
ordered.cols <- unlist(lapply(biodb, is.ordered))
biodb.ints <- biodb
biodb.ints[, ordered.cols] <- lapply(biodb.ints[, ordered.cols], function(x) as.numeric(x))

X.dummy.ints <- predict(dummyVars(prog.form, data = biodb.ints, fullRank = TRUE), newdata = biodb)

xgb.fit <- train(
    X.dummy,
    y,
    trControl = train.control,
    tuneGrid = tune.grid,
    metric = metric,
    method = "xgbTree"
)

xgb.fit.int <- train(
    X.dummy.ints,
    y,
    trControl = train.control,
    tuneGrid = tune.grid,
    metric = metric,
    method = "xgbTree"
)


gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                        n.trees = (1:30)*50,
                        shrinkage = c(0.001, 0.01, 0.1),
                        n.minobsinnode = 20)

gbm.fit <- train(
    X,
    y,
    method = 'gbm',
    metric = metric,
    trControl = train.control.up,
    tuneGrid = gbmGrid,
    verbose = FALSE
)

model.list <- list(GBM = gbm.fit, XGB = xgb.fit, XGB.int = xgb.fit.int, rpart = rpart.fit)
resamps <- resamples(model.list)
summary(resamps)

trellis.par.set(caretTheme())
dotplot(resamps, metric = c("F", "Sens", "Spec", "Accuracy"), main = "Progressed To Treatment Model Comparision")

model.varImp <- lapply(model.list, function(x) varImp(x)$importance)
all.varnames <- unique(unlist(lapply(model.varImp, rownames)))
var.importance <- data.frame(variable = all.varnames)
rownames(var.importance) <- all.varnames
for(model.name in names(model.varImp)) {
    imp.df <- model.varImp[[model.name]]
    var.importance[[model.name]] <- NA
    var.importance[rownames(imp.df), model.name] <- imp.df$Overall
}
var.importance

# do.call(cbind, lapply(model.varImp, `[[`, 'importance'))
confusionMatrix(predict(xgb.fit), reference = y, positive = 'yes')

plot(gbm.fit)
rpart.plot(rpart.fit$finalModel)
