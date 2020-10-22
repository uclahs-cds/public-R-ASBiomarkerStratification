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
    'BiopsyUpgraded', # This needs some work... do we want to predict the missing values as well?
    'Prostatectomy',
    'ProgressedToTreatment'
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
