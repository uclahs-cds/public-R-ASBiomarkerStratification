## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------
library(BoutrosLab.ASBiomarkerSynergy);
library(BoutrosLab.statistics.survival);
library(survival);
library(BoutrosLab.plotting.general);
library(flextable);
library(magrittr);
library(pROC);
library(caret);
library(here);
library(gtsummary);

knitr::opts_chunk$set(echo = FALSE, warning = FALSE, fig.align = 'center', out.width = '60%', fig.pos = 'H');

#axis.cex <- .7
#lab.cex <- .8
strip.cex <- 0.45

as.data <- default.load.data();
biodb <- as.data$biodb;

biomarkers <- load.biomarker.categories()

biomarkers.clinically.useful <- biomarkers[biomarkers$clinically.useful == 1,]

valid.patients <- !is.na(biodb$BiopsyUpgraded) #biodb$NoUpgradeAndProgressed == 0 | is.na(biodb$NoUpgradeAndProgressed);

biodb.valid <- biodb[valid.patients, ];

col.labels <- label.or.name(biodb)
names(col.labels) <- colnames(biodb)

seed <- 9999


## ----bx-biomarker-summary-------------------------------------------------------------------------------------------------------
demo.vars <- c(
            'Age',
            'BMI',
            'Race',
            'Ethnicity',
            'ProstateVolume',
            'SOCPSA',
            'PCA3',
            'T2ERG',
            'MiPSCancerRisk',
            'PercentFreePSA',
            'PHI',
            'GeneticRiskScore',
            'RSIlesionSignal',
            'RSIlesionPIRADS',
            'PSADensity',
            'PHIDensity'
            );

demo.vars <- c(
  'Age',
  'Race',
  'Ethnicity',
  'Weight',
  'Height',
  'BMI',
  'SOCPSA',
  'PSADensity',
  'freePSA', 'p2PSA', 'PercentFreePSA', 'PHI', 'PHIDensity',
  'MRIResult', 'MRILesions',
  'BiopsyResult',
  'ProstateVolume',
  'PreviousISUP',
  'StudyHighestISUP',
  'PCA3', 'T2ERG', 'MyProstateScore', 'MPSDensity',
  'GeneticRiskScore',
  'TNFaAverage',
  'RSInormalSignal',
  'RSIlesionSignal',
  'ADCnormalSignal',
  'ADClesionSignal',
  'RSIlesionPIRADS',
  'Germline.variants',
  'FollowUpTime'
)

target.vars <- c('ProgressedToTreatment', 'Prostatectomy')

demo.table <- tbl_summary(data = biodb[, c(demo.vars, 'BiopsyUpgraded')],
  by = BiopsyUpgraded,
  missing = 'no',
  label = RSIlesionPIRADS ~ 'RSI lesion PI-RADS') %>%
  add_p()

# demo.metadata <- demo.table$meta_data

target.table <- tbl_summary(data = biodb[, c(target.vars, 'BiopsyUpgraded')],
  by = BiopsyUpgraded,
  missing = 'no')

all.demo.table <- tbl_stack(list(demo.table, target.table)) %>%
  bold_labels() %>%
  modify_spanning_header(c('stat_1', 'stat_2') ~ '**Grade Group Upgrade**')

all.demo.table

flextable::save_as_docx(as_flex_table(all.demo.table), path = here::here('results/tables/demographics.docx'))


## -------------------------------------------------------------------------------------------------------------------------------
set.seed(seed)
create.heatmap.AS(biodb, filename = here('results/figures/corr_multiplot.tiff'), resolution = 300)


## -------------------------------------------------------------------------------------------------------------------------------
set.seed(seed)
uni.auc.ci <- lapply(biodb[, cor.variables], function(predictor) {
    pROC::ci.auc(
        roc(response = biodb$BiopsyUpgraded,
            predictor = predictor,
            plot = FALSE,
            direction = '<',
            levels = c('no', 'yes')))
})

count <- 0
for (i in names(uni.auc.ci)) {
  conf.int <- uni.auc.ci[[i]]
  if (conf.int[1] > 0.5) {
    count <- count + 1
    print(i)
    print(round(conf.int[c(2,1,3)], 2))
  }
}

targets <- c(
    'BiopsyUpgraded'#,
    #'Prostatectomy',
    #'ProgressedToTreatment'
    );
metric <- 'PR-AUC';
method <- 'gbm';
model.names <- c('everything', 'reduced_both', 'reduced_only_RSIlesionSignal', 'reduced_only_RSIlesionPIRADS');
model.nice.names <- c(
  'Full Model',
  'CR + RSI lesion signal and PI-RADS',
  'CR + RSI lesion signal',
  'CR + PI-RADS'
)
models <- lapply(targets, function(x) {
    model.id <- paste(x, metric, seed, sep = '_');
    model.id <- paste(model.id, model.names, sep = '_');
    model.files <- here(paste('models/gbm', model.id, 'model.RDS', sep = '_'));
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

optimal.thresholds <- lapply(models.roc, function(rocs) {
  lapply(rocs, function(r) {
      as.numeric(coords(r, 'best')[1])
  })
})


## -------------------------------------------------------------------------------------------------------------------------------
bx.preds <- rep(NA, nrow(biodb.valid))
bx.preds[!is.na(biodb.valid$BiopsyUpgraded)] <- ifelse(predict(models$BiopsyUpgraded$`Full Model`, type = 'prob')[, 'yes'] > optimal.thresholds$BiopsyUpgraded$`Full Model`, 1, 0)

bx.preds.factor <- as.factor(bx.preds)
levels(bx.preds.factor) <- c('no', 'yes')

bx.reduced.preds <- rep(NA, nrow(biodb.valid))
bx.reduced.preds[!is.na(biodb.valid$BiopsyUpgraded)] <- ifelse(predict(models$BiopsyUpgraded$`CR + RSI lesion signal`, type = 'prob')[, 'yes'] > optimal.thresholds$BiopsyUpgraded$`Full Model`, 1, 0)

bx.reduced.preds.factor <- as.factor(bx.reduced.preds)
levels(bx.reduced.preds.factor) <- c('no', 'yes')

confusionMatrix(bx.preds.factor, biodb.valid$BiopsyUpgraded)
confusionMatrix(bx.reduced.preds.factor, biodb.valid$BiopsyUpgraded)


## ----summary-table-bx-----------------------------------------------------------------------------------------------------------
summary.df <- summarize.models(models['BiopsyUpgraded'], models.roc['BiopsyUpgraded'])
cols <- c('model', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'F1')

thresholds.table <- summary.df[, cols]
# Break camel case into newlines
# thresholds.table$target <- kableExtra::linebreak(camel.to.spaces(thresholds.table$target, replace = '\n'), align = 'c')
num.cols <- 2:length(cols)
thresholds.table[, num.cols] <- lapply(thresholds.table[, num.cols], function(x) round(as.numeric(x), 2))

(upgrade.table <- flextable(thresholds.table) %>%
  add_footer_lines('CR = Clinically relevant variables') %>%
  # set_caption('Models comparison between the full model, and clinically relevant models with RSI lesion signal and PI-RADS. Including RSI lesion signal improves the model greatly while adding RSI PI-RADS does not have much effect if RSI lesion signal is in the model. All of the statistics are over the cross-validation 10-folds and 5 repetitions (N = 480).') %>%
  autofit())

save_as_docx(upgrade.table, path = here('results/tables/upgrade_table.docx'))


## -------------------------------------------------------------------------------------------------------------------------------
group.tests <- FALSE

if (group.tests) {
  seq.model.index <- c(1,4,5,6);
} else {
  seq.model.index <- 1:6
}

seq.model.id <- paste('sequential6', 'BiopsyUpgraded', metric, seed, sep = '_');
seq.model.file <- paste(seq.model.id, seq.model.index, 'model.RDS', sep = '_');
seq.model.path <- here::here(paste0('models/sequential/', seq.model.file));

seq.models <- lapply(seq.model.path, readRDS);

if (group.tests) {
  names(seq.models) <- c('Demographics', 'Blood, Urine, Genetics', 'MRI Features', 'Volume Corrected');
} else {
  names(seq.models) <- c('Demographics', 'Blood', 'Urine', 'Genetics', 'MRI Features', 'Volume Corrected');
}


seq.models.roc <- lapply(seq.models, function(m) {
    bestPreds <- with(m, merge(pred, bestTune));
    pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
})

best.thres <- lapply(seq.models.roc, function(x) {
  as.numeric(coords(x, 'best', transpose = TRUE)[1])
})
