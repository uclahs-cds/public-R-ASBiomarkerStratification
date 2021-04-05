library(BoutrosLab.ASBiomarkerSynergy);
# Baseline models
library(magrittr)
library(knitr)
library(kableExtra)
seed <- 1313;
targets <- c(
  'BiopsyUpgraded',
  'Prostatectomy',
  'ProgressedToTreatment');
metric <- 'PR-AUC';
method <- 'gbm';
model.names <- c('baseline', 'baseline_RSIlesionSignal', 'baseline_RSIlesionPIRADS');
model.nice.names <- c(
  'Baseline',
  'Baseline + RSI lesion Signal',
  'Baseline + RSI lesion PI-RADS'
)
models <- lapply(targets, function(x) {
  model.id <- paste(x, metric, seed, sep = '_');
  model.id <- paste(model.id, model.names, sep = '_');

  model.files <- here::here(paste('models/gbm', model.id, 'model.RDS', sep = '_'));
  names(model.files) <- model.nice.names;
  lapply(model.files, readRDS);
})
names(models) <- targets;

models.roc <- lapply(models, function(x) {
  lapply(x, function(m) {
    bestPreds <- with(m, merge(pred, bestTune));
    pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
  })
})

summary.df <- summarize.models(models, models.roc)

cols <- c('target', 'model', 'threshold', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1')

thresholds.table <- summary.df[, cols]
# Break camel case into newlines
thresholds.table$target <- kableExtra::linebreak(camel.to.spaces(thresholds.table$target, replace = '\n'), align = 'c')
num.cols <- 4:length(cols)
thresholds.table[, num.cols] <- lapply(thresholds.table[, num.cols], function(x) round(as.numeric(x), 2))

kable(thresholds.table,
      digits = 2,
      row.names = FALSE,
      caption = 'Cross-validation summary statistics for the different targets and models comparing RSI lesion signal and RSI lesion PI-RADS',
      escape = F) %>%
  kable_styling() %>%
  column_spec(1:2, width_min = '7em', width = '7em') %>%
  collapse_rows(columns = 1:2)

bx.upgrade.var.imp <- var.imp.combine(models$BiopsyUpgraded)

kable(x = bx.upgrade.var.imp[order(bx.upgrade.var.imp[,2], decreasing = TRUE), ],
      digits = 2,
      row.names = FALSE,
      caption = 'Variable importance for Biopsy Upgraded') %>%
  kable_styling() %>%
  column_spec(2:5, width_min = '5em', width = '5em')
