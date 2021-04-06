library(BoutrosLab.ASBiomarkerSynergy);
# Baseline models
library(magrittr)
library(knitr)
library(kableExtra)
seed <- 9999;
targets <- 'BiopsyUpgraded';
metric <- 'PR-AUC';
method <- 'gbm';
model.names <- c('baseline', 'baseline_RSIlesionSignal', 'baseline_RSIlesionPIRADS', 'baseline_both');
model.nice.names <- c(
  'Baseline',
  'Baseline + RSI',
  'Baseline + PI-RADS',
  'Baseline + RSI + PI-RADS'
)

set.seed(seed);
train.model <- FALSE
if(train.model) {
  biodb <- default.load.data(onlyBiodb = TRUE);

  # Baseline models
  lapply(targets, function(tg) {
    lapply(metric,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = AS.train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           baseline.model = TRUE,
           suffix = 'baseline'
    );
  });

  # Baseline + RSI lesion signal
  lapply(targets, function(tg) {
    lapply(metric,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = AS.train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           baseline.model = TRUE,
           include.vars = 'RSIlesionSignal',
           suffix = 'baseline_RSIlesionSignal'
    );
  });

  # Baseline + RSI PI-RADS
  lapply(targets, function(tg) {
    lapply(metric,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = AS.train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           baseline.model = TRUE,
           include.vars = 'RSIlesionPIRADS',
           suffix = 'baseline_RSIlesionPIRADS'
    );
  });

  # Baseline + RSI + PI-RADS
  lapply(targets, function(tg) {
    lapply(metric,
           AS.models,
           biodb = biodb,
           target = tg,
           train.control = AS.train.control,
           predict.missing = FALSE,
           seed = seed,
           models = c('gbm'),
           rm.NoUpgradeAndProgressed = FALSE,
           baseline.model = TRUE,
           include.vars = c('RSIlesionSignal', 'RSIlesionSignal'),
           suffix = 'baseline_both'
    );
  });
}

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

cols <- c('model', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1')

thresholds.table <- summary.df[, cols]
# Break camel case into newlines
# thresholds.table$target <- kableExtra::linebreak(camel.to.spaces(thresholds.table$target, replace = '\n'), align = 'c')
num.cols <- 2:length(cols)
thresholds.table[, num.cols] <- lapply(thresholds.table[, num.cols], function(x) round(as.numeric(x), 2))

(baseline.table <- flextable(thresholds.table) %>%
    add_footer_lines('Baseline = Age, BMI, Prostate Volume, PSA Density and %Free PSA. Models comparison between the full model, and clinically relevant models with RSI lesion signal and PI-RADS. Including RSI lesion signal improves the model greatly while adding RSI PI-RADS does not have much effect if RSI lesion signal is in the model. All of the statistics are over the cross-validation 10-folds and 5 repetitions. The optimal threshold is selected based on maximum sensitivity + specificity – 1).Baseline = Age, BMI, Prostate Volume, PSA Density and %Free PSA. Models comparison between the full model, and clinically relevant models with RSI lesion signal and PI-RADS. Including RSI lesion signal improves the model greatly while adding RSI PI-RADS does not have much effect if RSI lesion signal is in the model. All of the statistics are over the cross-validation 10-folds and 5 repetitions. The optimal threshold is selected based on maximum sensitivity + specificity – 1)') %>%
    set_caption('Predictive Models for Prostate Cancer Upgrading on Active Surveillance to compare RSI lesion signal versus RSI PI-RADS classification.') %>%
    autofit())

save_as_docx(baseline.table, path = here('results/tables/baseline_table.docx'))

