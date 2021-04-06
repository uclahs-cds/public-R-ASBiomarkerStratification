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
