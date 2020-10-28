library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

seed <- 1313;

train.control <- caret::trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5
    );

metrics <- c(#'ROC-AUC'
             'PR-AUC'
             );
targets <- c(
    'BiopsyUpgraded'#,
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
           rm.NoUpgradeAndProgressed = TRUE,
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
              rm.NoUpgradeAndProgressed = TRUE,
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
           rm.NoUpgradeAndProgressed = TRUE,
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
           rm.NoUpgradeAndProgressed = TRUE,
           reduced.model = FALSE,
           suffix = 'everything'
    );
});

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
#            rm.NoUpgradeAndProgressed = TRUE,
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
#            rm.NoUpgradeAndProgressed = TRUE,
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
#            rm.NoUpgradeAndProgressed = TRUE,
#            baseline.model = TRUE,
#            include.vars = 'RSIlesionPIRADS',
#            suffix = 'baseline_RSIlesionPIRADS'
#     );
# });
#
