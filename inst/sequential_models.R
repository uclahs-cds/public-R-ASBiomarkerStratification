library(BoutrosLab.ASBiomarkerSynergy);
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
# X <- lapply(seq(1, max(useful.biomarkers$order)), function(i) {
#   bio.vars <- useful.biomarkers$variable[useful.biomarkers$order <= i]
#   biodb[!missing.target & valid.patients, bio.vars, drop = FALSE]
# })

# Separate blood, urine, genetics
biocategories <- unique(biomarkers$category)

X <- lapply(seq_along(biocategories), function(i) {
  bio.cats <- biocategories[1:i];
  bio.vars <- useful.biomarkers$variable[useful.biomarkers$category %in% bio.cats]
  biodb[!missing.target & valid.patients, bio.vars, drop = FALSE]
})

y.target <- biodb[!missing.target & valid.patients, target];
y <- y.target
levels(y) <- c('no', 'yes');

# Save the models to file
# seq.id <- c('demographics', 'pre-MRI', 'MRI', 'Post-MRI');
if(length(X) == 6) {
  model.id <- paste('sequential6', target, metric, seed, sep = '_');
} else {
  model.id <- paste('sequential4', target, metric, seed, sep = '_');
}

model.file <- paste(model.id, 1:length(X), 'model.RDS', sep = '_');

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

if(length(gbm.seq.models) == 4) {
  names(gbm.seq.models) <- c('Demographics', 'Blood/Urine/Genetics', 'MRI Features', 'Post MRI');
}
if(length(gbm.seq.models) == 6) {
  names(gbm.seq.models) <- c('Demographics', 'Blood', 'Urine', 'Genetics', 'MRI Features', 'Post MRI');
}



