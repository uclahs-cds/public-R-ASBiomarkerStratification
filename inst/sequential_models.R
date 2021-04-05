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
seed <- 9999
missing.target <- is.na(biodb[, target]);

biomarkers.clinically.useful <- biomarkers[biomarkers$clinically.useful == 1, ];

# Separate blood, urine, genetics
biocategories <- unique(biomarkers$category)

col.labels <- label.or.name(biodb)
names(col.labels) <- colnames(biodb)

X <- lapply(seq_along(biocategories), function(i) {
  bio.cats <- biocategories[1:i];
  print(bio.cats);
  bio.vars <- biomarkers.clinically.useful$variable[biomarkers.clinically.useful$category %in% bio.cats]
  print(bio.vars);
  biodb[!missing.target, bio.vars, drop = FALSE]
})

y.target <- biodb[!missing.target, target];
y <- y.target
levels(y) <- c('no', 'yes');

# Save the models to file
# seq.id <- c('demographics', 'pre-MRI', 'MRI', 'Post-MRI');
if (length(X) == 6) {
  model.id <- paste('sequential6', target, metric, seed, sep = '_');
} else {
  model.id <- paste('sequential4', target, metric, seed, sep = '_');
}

model.file <- paste(model.id, 1:length(X), 'model.RDS', sep = '_');

seq.models <- lapply(seq_along(X), function(i) {
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

if (length(seq.models) == 4) {
  names(seq.models) <- c('Demographics', 'Blood/Urine/Genetics', 'MRI Features', 'Volume Corrected');
}
if (length(seq.models) == 6) {
  names(seq.models) <- c('Demographics', 'Blood', 'Urine', 'Genetics', 'MRI Features', 'Volume Corrected');
}

# Load instead of training
seq.model.id <- paste('sequential6', 'BiopsyUpgraded', metric, seed, sep = '_');
seq.model.file <- paste(seq.model.id, 1:6, 'model.RDS', sep = '_');
seq.model.path <- here::here(paste0('models/sequential/', seq.model.file));

seq.models <- lapply(seq.model.path, readRDS);

seq.models.roc <- lapply(seq.models, function(m) {
  bestPreds <- with(m, merge(pred, bestTune));
  pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
})

best.thres <- lapply(seq.models.roc, function(x) {
  as.numeric(coords(x, 'best', transpose = TRUE)[1])
})

seq.var.names <- lapply(seq.models, function(x) x$finalModel$var.names)
all.var.names <- biomarkers.clinically.useful$variable

summary.df <- summarize.seq.models(seq.models, seq.models.roc)
metrics <- c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1")
cols <- c("group", metrics)

thresholds.table <- summary.df[, cols]
thresholds.table$group <- lapply(seq_along(seq.models), function(i) {
  paste0(names(seq.models)[1:i], collapse = " +\n")
})

# Create the dotmap
dotmap.data <- matrix(0, ncol = length(seq.models), nrow = length(all.var.names))
colnames(dotmap.data) <- 1:length(seq.models)
rownames(dotmap.data) <- all.var.names

for (i in seq_along(seq.models)) {
  dotmap.data[seq.var.names[[i]], i] <- 1
}

rownames(dotmap.data) <- col.labels[all.var.names]

var.categories <- factor(biomarkers.clinically.useful$category, levels = unique(biomarkers.clinically.useful$category))

dotmap.bg.data <- as.data.frame(matrix(rep(var.categories, 6), nrow = length(biomarkers.clinically.useful$category), ncol = 6))

dotmap.bg.data <- lapply(dotmap.bg.data, factor, levels = unique(biomarkers.clinically.useful$category))
dotmap.bg.data.int <- data.frame(lapply(dotmap.bg.data, as.numeric))

colnames(dotmap.bg.data.int) <- 1:length(seq.models)
rownames(dotmap.bg.data.int) <- all.var.names

bpg.dotmap <- create.dotmap(
  x = dotmap.data,
  spot.size.function = identity,
  spot.colour.function = function(x) "black",
  colour.scheme = default.colours(nlevels(var.categories), palette.type = 'qual'),
  total.colours = nlevels(var.categories) + 1,
  bg.alpha = 0.65,
  bg.data = dotmap.bg.data.int
)

metric.scatterplot.data <- thresholds.table %>%
  dplyr::mutate(group = 1:length(seq.models)) %>%
  tidyr::gather(metric, value, -group)
#
metric.colours <- c('dodgerblue', 'goldenrod1', 'darkorange1', 'seagreen2', 'orchid3')
names(metric.colours) <- metrics

metric.scatterplot <- create.scatterplot(
  formula = value ~ group,
  data = metric.scatterplot.data,
  groups = metric.scatterplot.data$metric,
  left.padding = 0,
  xat = seq_along(seq.models),
  col = metric.colours,
  type = 'o',
  add.grid = TRUE,
  # Remove the vertial grid lines
  xgrid.at = 0,
  cex = 2,
  lwd = 4
)

metric.legend <- legend.grob(
  list(
    legend = list(
      colours = metric.colours,
      title = "Metrics",
      labels = metrics,
      border = 'black'
    ),
    legend = list(
      colours = default.colours(nlevels(var.categories), palette.type = 'qual'),
      title = "Test Methodology",
      labels = levels(var.categories),
      border = 'black'
    )
  ),
  title.cex = 2,
  label.cex = 2,
  title.just = 'left'
);

create.multiplot(
  plot.objects = list(bpg.dotmap, metric.scatterplot),
  panel.heights = c(1, 1),
  y.relation = 'free',
  left.padding = 30,
  right.padding = 2,
  print.new.legend = TRUE,
  xat = seq_along(seq.models),
  xaxis.labels = seq_along(seq.models),
  xlab.label = 'Sequence Group',

  ylab.label = list(
    'Metric Value',
    ''
  ),
  legend = list(
    right = list(fun = metric.legend)),
  height = 12,
  width = 16,
  resolution = 300,
  filename = here('euro_urology/figures/seq_dotmap.tiff')
)
