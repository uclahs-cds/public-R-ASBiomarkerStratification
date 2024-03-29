# To determine the optimal combination of these models, we first formed a model of clinico-epidemiologic variables only.
# This model has an AUROC of XX (95% CI: X-Y) and an AUPRC of YY (95% CI: Z-Q; Figure 1D).
# We then generated three biomarker-driven models:
#  - one with clinico-epidemiologic and radiologic features
#  - one with clinico-epidemiologic and molecular features
#  - one containing all three

library(BoutrosLab.ASBiomarkerSynergy);
library(caret);
library(pROC);

seed <- 9999;
metric <- 'PR-AUC';
target <- 'BiopsyUpgraded';

as.data <- default.load.data();
biodb <- as.data$biodb;
biomarkers <- load.biomarker.categories();
rownames(biomarkers) <- biomarkers$variable;
biomarkers <- biomarkers[cor.variables,]

# baseline.vars <- biomarkers$variable[biomarkers$category == 'Demographics']
clinico.epi.vars <- c('Age', 'BMI', 'ProstateVolume', 'PSADensity', 'PercentFreePSA');
radiologic.vars <- c(clinico.epi.vars, biomarkers$variable[biomarkers$category == 'MRI Features']);
molecular.vars <- c(clinico.epi.vars, biomarkers$variable[biomarkers$category == 'Genetics']);
all.vars <- union(radiologic.vars, molecular.vars);

train.model <- FALSE;
if (train.model) {
  dir.create(here::here('models'), showWarnings = FALSE, recursive = TRUE);
  missing.target <- is.na(biodb[, target]);
  X.clinico.epi <- biodb[!missing.target, clinico.epi.vars]
  X.radiologic <- biodb[!missing.target, radiologic.vars]
  X.molecular <- biodb[!missing.target, molecular.vars]
  X <- biodb[!missing.target, all.vars];

  y <- biodb[!missing.target, target];

  set.seed(seed);

  gbm.grid <- gbm.hyper.grid();

  print('Fitting gbm model for clinico-epidemiologic features');
  gbm.clinico.epi <- caret::train(
    X.clinico.epi,
    y,
    method = 'gbm',
    metric = metric,
    trControl = AS.train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
  )
  print('Completed fitting gbm model...');

  model.id <- paste(target, metric, seed, sep = '_');
  model.file <- paste('gbm', model.id, 'clinico_epidemiologic_model.RDS', sep = '_');
  # print(paste0('Saving file to: ',  here::here(paste0('models/', model.file))));
  saveRDS(object = gbm.clinico.epi, file = here::here(paste0('models/', model.file)));


  print('Fitting gbm model for radiologic  features');
  gbm.radiologic <- caret::train(
    X.radiologic,
    y,
    method = 'gbm',
    metric = metric,
    trControl = AS.train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
  )
  print('Completed fitting gbm model...');

  model.id <- paste(target, metric, seed, sep = '_');
  model.file <- paste('gbm', model.id, 'radiologic_model.RDS', sep = '_');
  # print(paste0('Saving file to: ',  here::here(paste0('models/', model.file))));
  saveRDS(object = gbm.radiologic, file = here::here(paste0('models/', model.file)));

  print('Fitting gbm model for molecular features');
  gbm.molecular <- caret::train(
    X.molecular,
    y,
    method = 'gbm',
    metric = metric,
    trControl = AS.train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
  )
  print('Completed fitting gbm model...');

  model.id <- paste(target, metric, seed, sep = '_');
  model.file <- paste('gbm', model.id, 'molecular_model.RDS', sep = '_');
  saveRDS(object = gbm.molecular, file = here::here(paste0('models/', model.file)));

  print('Fitting gbm model for all features');
  gbm.all <- caret::train(
    X,
    y,
    method = 'gbm',
    metric = metric,
    trControl = train.control,
    tuneGrid = gbm.grid,
    verbose = FALSE
    );
  print('Completed fitting gbm model...');

  model.id <- paste(target, metric, seed, sep = '_');
  model.file <- paste('gbm', model.id, 'combine_model.RDS', sep = '_');
  saveRDS(object = gbm.all, file = here::here(paste0('models/', model.file)));
  } else {
  model.id <- paste(target, metric, seed, sep = '_');

  # Clinico
  model.file <- paste('gbm', model.id, 'clinico_epidemiologic_model.RDS', sep = '_');
  gbm.clinico.epi <- readRDS(file = here::here(paste0('models/', model.file)));

  # radiologic
  model.file <- paste('gbm', model.id, 'radiologic_model.RDS', sep = '_');
  gbm.radiologic <- readRDS(file = here::here(paste0('models/', model.file)));

  # molecular
  model.file <- paste('gbm', model.id, 'molecular_model.RDS', sep = '_');
  gbm.molecular <- readRDS(file = here::here(paste0('models/', model.file)));

  # All variables model
  model.file <- paste('gbm', model.id, 'combine_model.RDS', sep = '_');
  gbm.all <- readRDS(file = here::here(paste0('models/', model.file)));
  }


set.seed(seed);
models <- list(gbm.clinico.epi, gbm.radiologic, gbm.molecular, gbm.all);
names(models) <- c('clinico-epidemiologic', 'radiologic', 'molecular', 'all');

models.roc <- lapply(models, function(m) {
  bestPreds <- with(m, merge(pred, bestTune));
  pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
  });

models.cvRoc <- lapply(models, function(m) {
  bestPreds <- with(m, merge(pred, bestTune));
  cvAUC::ci.cvAUC(predictions = bestPreds$yes, labels = bestPreds$obs, folds = bestPreds$Resample);
  });

auc.ci.format <- lapply(models.cvRoc, function(r) {
  sprintf('%.2f (95%% CI: %.2f - %.2f)', r$cvAUC, r$ci[1], r$ci[2]);
  });

auc.ci.format

models.pr <- lapply(models, function(m) {
  bestResults <- with(m, merge(results, bestTune));
  round(bestResults[c('PR-AUC', 'PR-AUCSD')], 2)
})
models.pr

# Save the ROC/PR Curves
tiff(filename = here::here('results/figures/roc-pr_curves_alt.tiff'), width = 10, height = 8, res = 300, units = 'in')
roc.pr.plot(models.roc)
dev.off()

models.target <- list('BiopsyUpgraded' = models);
models.roc.target <- list('BiopsyUpgraded' = models.roc);
summary.df <- summarize.models(models.target, models.roc.target);
cols <- c('model', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1');

thresholds.table <- summary.df[, cols];
num.cols <- 2:length(cols);
thresholds.table[, num.cols] <- lapply(thresholds.table[, num.cols], function(x) round(as.numeric(x), 2));

# Print and save the results table
(upgrade.table <- flextable::autofit(flextable::flextable(thresholds.table)));
flextable::save_as_docx(upgrade.table, path = here::here('results/tables/biomarker_table.docx'))
