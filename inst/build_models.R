library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

train.control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 10
    );

metrics <- c('F');
target <- c(
#    'ProgressedToTreatment',
     'BiopsyUpgraded', # This needs some work... do we want to predict the missing values as well?
    'Prostatectomy');

results <- lapply(target, eval.models, biodb = biodb, train.control = train.control, metrics = metrics, predict.missing = FALSE, seed = 1313);


