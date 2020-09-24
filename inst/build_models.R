library(ProstateCancer.ASBiomarkerSynergy);

biodb <- default.load.data(onlyBiodb = TRUE);

train.control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 10
    );

metrics <- c('F', 'Accuracy');
target <- c('ProgressedToTreatment', 'BiopsyUpgraded', 'Prostatectomy');

results <- lapply(target, eval.models, biodb = biodb, train.control = train.control, metrics = metrics);
length(results)
results[[1]]$resamples

models <- results[[1]]$models

model.resamples <- lapply(seq_along(metrics), function(i) {
    metric <- metrics[i]
    model <- models[[i]]
    model.list <- list(model$gbm.fit, model$c50.fit, model$rpart.fit)
    names(model.list) <- paste0(c('GBM', 'C5.0', 'rpart'), '.', metric)
    resamples(model.list)
})

results[[1]]$resamples <- model.resamples

models <- results[[3]]$models

model.resamples <- lapply(seq_along(metrics), function(i) {
    metric <- metrics[i]
    model <- models[[i]]
    model.list <- list(model$gbm.fit, model$c50.fit, model$rpart.fit)
    names(model.list) <- paste0(c('GBM', 'C5.0', 'rpart'), '.', metric)
    resamples(model.list)
})

results[[3]]$resamples <- model.resamples
