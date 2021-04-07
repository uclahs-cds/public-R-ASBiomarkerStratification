library(BoutrosLab.ASBiomarkerSynergy);

seed <- 9999;
metric <- 'PR-AUC'
target <- 'BiopsyUpgraded'
train.model <- FALSE

as.data <- default.load.data();
biodb <- as.data$biodb;
biomarkers <- load.biomarker.categories();
rownames(biomarkers) <- biomarkers$variable;
biomarkers <- biomarkers[cor.variables,]

biodb.valid <- biodb[!is.na(biodb$BiopsyUpgraded), ];

clinico.epi.vars <- c('Age', 'BMI', 'ProstateVolume', 'PSADensity', 'PercentFreePSA');
radiologic.vars <- c(clinico.epi.vars, biomarkers$variable[biomarkers$category == 'MRI Features']);
molecular.vars <- c(clinico.epi.vars, biomarkers$variable[biomarkers$category == 'Genetics'])
all.vars <- union(radiologic.vars, molecular.vars);

model.id <- paste(target, metric, seed, sep = '_');

# Clinico
model.file <- paste('gbm', model.id, 'clinico_epidemiologic_model.RDS', sep = '_');
gbm.clinico.epi <- readRDS(file = here::here(paste0('models/', model.file)));

# All variables model
model.file <- paste('gbm', model.id, 'combine_model.RDS', sep = '_');
gbm.all <- readRDS(file = here::here(paste0('models/', model.file)));

set.seed(seed);
models <- list(gbm.clinico.epi, gbm.all);
names(models) <- c('clinico-epidemiologic', 'all');

models.roc <- lapply(models, function(m) {
  bestPreds <- with(m, merge(pred, bestTune));
  pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
  });

optimal.thresholds <- lapply(models.roc, function(r) {
    as.numeric(coords(r, 'best', transpose = FALSE)[1])
  });

bx.preds <- rep(NA, nrow(biodb.valid));
bx.num.preds <- predict(models$BiopsyUpgraded$`Full Model`, type = 'prob')[, 'yes'];
bx.preds[!is.na(biodb.valid$BiopsyUpgraded)] <- ifelse(bx.num.preds > optimal.thresholds$BiopsyUpgraded$`Full Model`, 1, 0);

bx.preds.factor <- as.factor(bx.preds);
levels(bx.preds.factor) <- c('no', 'yes');

bx.clinico.epi.preds <- ifelse(predict(models$`clinico-epidemiologic`, type = 'prob')[, 'yes'] > optimal.thresholds$`clinico-epidemiologic`, 1, 0);
bx.clinico.epi.preds.factor <- as.factor(bx.clinico.epi.preds);
levels(bx.clinico.epi.preds.factor) <- c('no', 'yes');

bx.all.preds <- ifelse(predict(models$all, type = 'prob')[, 'yes'] > optimal.thresholds$all, 1, 0);
bx.all.preds.factor <- as.factor(bx.all.preds);
levels(bx.all.preds.factor) <- c('no', 'yes');

confusionMatrix(bx.clinico.epi.preds.factor, biodb.valid$BiopsyUpgraded, positive = 'yes');
confusionMatrix(bx.all.preds.factor, biodb.valid$BiopsyUpgraded, positive = 'yes');

axis.cex <- 1.1;
lab.cex <- 1.6;

km.arguments <- list(
  xaxis.cex = axis.cex,
  yaxis.cex = axis.cex,
  xlab.cex = lab.cex,
  ylab.cex = lab.cex,
  show.risktable = TRUE,
  risktable.fontsize = 10,
  key.groups.cex = 1.3,
  key.stats.cex = 1.3,
  key.stats.y.pos = 0.5,
  return.statistics = FALSE,
  digits = 1
  );

upgrade.surv.valid <- do.call(
  what = survival::Surv,
  args = surv.format(biodb.valid$DaysDxToUpgrade, biodb.valid$DaysDxToLastReview)
  );

patient.groups <- as.factor(bx.all.preds);
levels(patient.groups) <- c('No upgrade prediction', 'Predict upgrade');

survival::survfit(upgrade.surv.valid ~ patient.groups);
summary(survival::survfit(upgrade.surv.valid ~ patient.groups));

# Create the survival plot
do.call(
  what = BoutrosLab.plotting.survival::create.km.plot,
  args = c(
    km.arguments,
    list(
      main.cex = 1.75,
      survival.object = upgrade.surv.valid,
      patient.groups = patient.groups,
      main = 'Overall Days-to-Upgrade from Diagnosis',
      xlab.label = 'Days to upgrade since diagnosis',
      ylab.label = 'Upgrade-free survival',
      filename = here('results/figures/km_dx-to-upgrade_v2.tiff'),
      left.padding = 10,
      risk.label.pos = -1600,
      resolution = 1000
      )
    )
  );


# Clinico-Epi TTE
patient.groups.epi <- bx.clinico.epi.preds.factor;
levels(patient.groups.epi) <- c('No upgrade prediction', 'Predict upgrade');

survival::survfit(upgrade.surv.valid ~ patient.groups.epi);
summary(survival::survfit(upgrade.surv.valid ~ patient.groups.epi));

do.call(
  what = BoutrosLab.plotting.survival::create.km.plot,
  args = c(
    km.arguments,
    list(
      main.cex = 1.75,
      survival.object = upgrade.surv.valid,
      patient.groups = patient.groups.epi,
      main = 'Overall Days-to-Upgrade from Diagnosis',
      xlab.label = 'Days to upgrade since diagnosis',
      ylab.label = 'Upgrade-free survival',
      filename = here('results/figures/km_dx-to-upgrade_epi_v2.tiff'),
      left.padding = 10,
      risk.label.pos = -1600
      )
    )
  );
