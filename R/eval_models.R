#' Title
#'
#' @param biodb
#' @param target
#' @param metrics
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
eval.models <- function(
    biodb,
    target = c('ProgressedToTreatment', 'BiopsyUpgraded', 'Prostatectomy'),
    metrics = c('Accuracy', 'F', 'AUC', 'Kappa', 'Precision', 'Recall', 'ROC',
                                           'Sens', 'Spec'),
    ...) {
    if(missing(metrics)) {
        metrics <- c('F', 'ROC');
        } else {
        metrics <- match.arg(metrics, several.ok = TRUE);
        }

    target <- match.arg(target);

    models <- lapply(metrics, AS.models, biodb = biodb, target = target, ...);

    # model.file <- paste(target, metrics, 'model.RDS', sep = "_")
    # # Save RDS of the models
    # lapply(seq_along(models), function(i) {
    #     saveRDS(object = models[[i]], file = here(paste0('models/', model.file[i])))
    # });

    # model.resamples <- lapply(seq_along(metrics), function(i) {
    #     metric <- metrics[i]
    #     model <- models[[i]]
    #     model.list <- list(model$gbm.fit, model$xgb.grid, model$rpart.fit)
    #     names(model.list) <- paste0(c('GBM', 'XGB', 'rpart'), '.', metric)
    #     resamples(model.list)
    # })
    #
    # list(
    #     models = models,
    #     resamples = model.resamples
    # )
}
