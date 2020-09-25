library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

train.control <- trainControl(
    method = "repeatedcv",
    number = 2,
    repeats = 2
    );

metrics <- c('F');
targets <- c(
    'ProgressedToTreatment',
    'BiopsyUpgraded', # This needs some work... do we want to predict the missing values as well?
    'Prostatectomy');

targets <- 'ProgressedToTreatment'

results <- lapply(targets, function(tg) {
    res <- lapply(metrics, AS.models, biodb = biodb, target = tg, train.control = train.control, predict.missing = FALSE, seed = 10101);
    names(res) <- metrics;
    res;
});
names(results) <- targets;

compare.var.imp(results$ProgressedToTreatment$F)
