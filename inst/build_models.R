library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

train.control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5
    );

metrics <- c('F');
targets <- c(
    'BiopsyUpgraded', # This needs some work... do we want to predict the missing values as well?
    'Prostatectomy',
    'ProgressedToTreatment');

results <- lapply(targets, function(tg) {
    res <- lapply(metrics, AS.models, biodb = biodb, target = tg, train.control = train.control, predict.missing = FALSE, seed = 1010);
    names(res) <- metrics;
    res;
});
names(results) <- targets;

compare.var.imp(results$ProgressedToTreatment$F, include.ranks = TRUE)
