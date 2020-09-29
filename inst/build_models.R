library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

seed <- 456456;

train.control <- caret::trainControl(
    method = "repeatedcv",
    number = 20,
    repeats = 5,
    sampling = "up"
    );

metrics <- c('PR-AUC');
targets <- c(
    'BiopsyUpgraded', # This needs some work... do we want to predict the missing values as well?
    'Prostatectomy',
    'ProgressedToTreatment');

results <- lapply(targets, function(tg) {
    res <- lapply(metrics,
                  AS.models,
                  biodb = biodb,
                  target = tg,
                  train.control = train.control,
                  predict.missing = FALSE,
                  seed = seed,
                  rpart.cost = matrix(c(0,1,2,0), byrow = TRUE, nrow = 2));
    names(res) <- metrics;
    res;
});
names(results) <- targets;

# compare.var.imp(results$ProgressedToTreatment$F, include.ranks = TRUE)
