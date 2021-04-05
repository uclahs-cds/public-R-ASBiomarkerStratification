library(BoutrosLab.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

seed <- 9999;

train.control <- caret::trainControl(
    method = 'repeatedcv',
    number = 10,
    repeats = 5
    );

metrics <- c(#'ROC-AUC'
             'PR-AUC'
             );
targets <- c(
    'BiopsyUpgraded' #,
    #'Prostatectomy',
    # 'ProgressedToTreatment'
    );

lapply(targets, function(tg) {
    lapply(metrics,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           reduced.model = TRUE,
           suffix = 'reduced_both'
    );
});

lapply(targets, function(tg) {
    lapply(metrics,
              AS.models,
              biodb = biodb,
              target = tg,
              train.control = train.control,
              predict.missing = FALSE,
              seed = seed,
              models = c('gbm'),
              rm.NoUpgradeAndProgressed = FALSE,
              reduced.model = TRUE,
              exclude.vars = 'RSIlesionSignal',
              suffix = 'reduced_only_RSIlesionPIRADS'
              );
    });


lapply(targets, function(tg) {
    lapply(metrics,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           reduced.model = TRUE,
           exclude.vars = 'RSIlesionPIRADS',
           suffix = 'reduced_only_RSIlesionSignal'
    );
});


lapply(targets, function(tg) {
    lapply(metrics,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           reduced.model = FALSE,
           suffix = 'everything'
    );
});

## Checking if adding additional variables to RSI lesion signal / PIRADS makes a difference
target <- 'BiopsyUpgraded'
demo.vars <- biomarkers[biomarkers$category == 'Demographics' & biomarkers$clinically.useful == 1, 'variable']
only.rsi.vars <- c(demo.vars, 'RSIlesionPIRADS', 'RSIlesionSignal')

missing.target <- is.na(biodb[, target]);
y.target <- biodb[!missing.target, target];
y <- y.target
levels(y) <- c('no', 'yes');
X <- biodb[!missing.target, only.rsi.vars];

gbm.grid <- gbm.hyper.grid()

train.control <- trainControl(
    method = 'repeatedcv',
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    savePredictions = T,
    summaryFunction = custom.summary
)

set.seed(seed);

print('Fitting gbm model');
gbm.rsi.fit <- caret::train(
    X,
    y,
    method = 'gbm',
    metric = metric,
    trControl = train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
)



# Drop the missing target values


# Baseline models
# lapply(targets, function(tg) {
#     lapply(metrics,
#            AS.models,
#            biodb = biodb,
#            target = tg,
#            train.control = train.control,
#            predict.missing = FALSE,
#            seed = seed,
#            models = c('gbm'),
#            rm.NoUpgradeAndProgressed = FALSE,
#            baseline.model = TRUE,
#            suffix = 'baseline'
#     );
# });
#
# # Baseline + RSI lesion signal
# lapply(targets, function(tg) {
#     lapply(metrics,
#            AS.models,
#            biodb = biodb,
#            target = tg,
#            train.control = train.control,
#            predict.missing = FALSE,
#            seed = seed,
#            models = c('gbm'),
#            rm.NoUpgradeAndProgressed = FALSE,
#            baseline.model = TRUE,
#            include.vars = 'RSIlesionSignal',
#            suffix = 'baseline_RSIlesionSignal'
#     );
# });
#
# # Baseline + RSI PI-RADS
# lapply(targets, function(tg) {
#     lapply(metrics,
#            AS.models,
#            biodb = biodb,
#            target = tg,
#            train.control = train.control,
#            predict.missing = FALSE,
#            seed = seed,
#            models = c('gbm'),
#            rm.NoUpgradeAndProgressed = FALSE,
#            baseline.model = TRUE,
#            include.vars = 'RSIlesionPIRADS',
#            suffix = 'baseline_RSIlesionPIRADS'
#     );
# });
#
