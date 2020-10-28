library(ProstateCancer.ASBiomarkerSynergy);
library(caret)
library(pROC)
library(rpart.plot)

train.control <- trainControl(
  method = 'repeatedcv',
  number = 10,
  repeats = 5,
  classProbs = TRUE,
  savePredictions = T,
  summaryFunction = custom.summary
)

gbm.grid <- gbm.hyper.grid()

biodb <- default.load.data(onlyBiodb = TRUE);
# roc.res <- roc(biodb$BiopsyUpgraded, predictor = biodb$RSIlesionSignal, algorithm = 0)

# Biomarkers and categories
biomarkers <- load.biomarker.categories()

target <- 'BiopsyUpgraded'
metric <- 'PR-AUC'
seed <- 1313
missing.target <- is.na(biodb[, target]);

# Only keep the patient that did not leave AS voluntarily if we have the rm.NoUpgradeAndProgressed flag
valid.patients <- (! (biodb$NoUpgradeAndProgressed == 1) | is.na(biodb$NoUpgradeAndProgressed));

useful.biomarkers <- biomarkers[biomarkers$clinically.useful == 1, ]

# List of data frame with the sequential data frames
X <- lapply(seq(1, max(useful.biomarkers$order)), function(i) {
  bio.vars <- useful.biomarkers$variable[useful.biomarkers$order <= i]
  biodb[!missing.target & valid.patients, bio.vars, drop = FALSE]
})

y.target <- biodb[!missing.target & valid.patients, target];
y <- y.target
levels(y) <- c('no', 'yes');

# Save the models to file
# seq.id <- c('demographics', 'pre-MRI', 'MRI', 'Post-MRI');
model.id <- paste('sequential', target, metric, seed, sep = '_');
model.file <- paste(model.id, 1:4, 'model.RDS', sep = '_');

gbm.seq.models <- lapply(seq_along(X), function(i) {
  x <- X[[i]];
  save.file <- model.file[[i]];
  set.seed(seed);
  mod <- caret::train(
    x,
    y,
    method = 'gbm',
    metric = 'PR-AUC',
    trControl = train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
  );

  saveRDS(object = mod, file = here::here(paste0('models/sequential/', save.file)));
  mod
})

names(gbm.seq.models) <- c('Demographics', 'Blood/Urine/Genetics', 'MRI Features', 'Post MRI + Biomarkers');

