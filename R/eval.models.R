#' Computes the optimal threshold based on a given cost matrix.
#'
#' @param model
#' @param cost.matrix
#' @param thresholds
#' @param smooth
#'
#' @return
#' @export
optimal.threshold.train <- function(model,
                                    cost.matrix = matrix(c(0,1,1,0), byrow = TRUE, nrow = 2),
                                    thresholds = c(0.25, 0.5, 0.75),
                                    smooth = TRUE) {
    threshold.cost <- cost.threshold.train(model, thresholds = thresholds, cost.matrix = cost.matrix)
    if (smooth) {
        threshold.cost <- as.vector(moving.avg(threshold.cost, order = max(1, round(thresholds / 20))))
    }
    thresholds[which.min(threshold.cost)]
}

flatten.ConfusionMatrix <- function(...) {
    res <- caret::confusionMatrix(...)
    cbind.data.frame(t(res$overall), t(res$byClass))
}

#' Compute the threshold summary stats per fold
#'
#' @param model
#' @param threshold
#'
#' @return
#' @export
threshold.summary.stats <- function(model, threshold, per.fold = FALSE) {
    pred.best <- with(model,
                      merge(pred, bestTune))
    pred.best$thres.preds <- as.factor(ifelse(pred.best[, 'yes'] > threshold, 'yes', 'no'))
    if (per.fold) {
        fold.split <- split(pred.best, pred.best$Resample)
        res <- lapply(fold.split, function(x) flatten.ConfusionMatrix(x$thres.preds, reference = x$obs, positive = 'yes', mode = 'everything'))
        do.call(rbind.data.frame, res)
    } else {
        flatten.ConfusionMatrix(pred.best$thres.preds, reference = pred.best$obs, positive = 'yes', mode = 'everything')
    }
}

# Combines the importance score of dummy variables into one
squash.importance <- function(importance, FUN = 'sum') {
    df <- importance
    df$row.names <- row.names(df)
    df$row.names <- gsub('(.*)\\..*', '\\1', df$row.names)

    agg.df <- aggregate(df$Overall, by = list(vars = df$row.names), FUN = FUN)
    res <- data.frame(Overall = agg.df[, 2])
    rownames(res) <- agg.df[, 1]
    res
}

#' Compute the variable importance for GBM, XGB, rpart
#'
#' @param x
#' @param onlyImportance boolean indicating whether only a data frame with the importance scores be returned? FALSE returns class 'caret::varImp.train'.
#' @param squash Should multiple dummy variables be squashed into one variable? (Relevant for xgb models)
#'
#' @return
#' @export
var.imp <- function(x, onlyImportance = TRUE, squash = TRUE) {
    stopifnot('train' == class(x));
    if ('gbm' == x$method) {
        gbm.sum <- gbm::summary.gbm(x$finalModel, plotit = FALSE);
        importance <- data.frame(Overall = gbm.sum[,2], row.names = gbm.sum[,1]);
        res <- list(
            model = 'gbm',
            calledFrom = 'summary.gbm',
            importance = importance
        );
        class(res) <- 'varImp.train';
    } else {
        res <- caret::varImp(x);
    }
    row.names(res$importance) <- gsub('(.*)\\.1', '\\1', row.names(res$importance));
    if (squash) {
        res$importance <- squash.importance(res$importance);
    }
    res$importance$variable <- row.names(res$importance);
    if (onlyImportance) res <- res$importance;
    res
}

#' Combine the variable importance scores from multiple models into one data frame
#'
#' @param x list of caret train classes
#'
#' @return
#' @export
var.imp.combine <-  function(x, order = FALSE, rev.cols = FALSE) {
    if (rev.cols) {
        x <- x[rev(names(x))];
        }
    res <- Reduce(function(...) merge(..., by = 'variable', all = TRUE), lapply(x, var.imp));
    rownames(res) <- res$variable;

    if (!is.null(names(x))) {
        colnames(res) <- c('variable', names(x));
    }

    if (order) {
        order.args <- as.list(res[, names(x)])
        order.args$decreasing <- TRUE
        res.order <- do.call('order', order.args)
        res <- res[res.order, ]
    }

    res
}

#' Compare variable importance across multiple models
#'
#' @param x
#' @param ranks
#'
#' @return
#' @export
compare.var.imp <- function(x, include.ranks = FALSE) {
    model.varImp <- lapply(x, function(x) var.imp(x)$importance)
    if (is.null(names(x))) {
        names(x) <- unlist(lapply(x, `[[`, 'method'));
    }

    all.varnames <- unlist(lapply(model.varImp, rownames));
    all.varnames <- unique(gsub('(.*)\\.1', '\\1', all.varnames));
    var.importance <- data.frame(variable = all.varnames);
    rownames(var.importance) <- all.varnames;
    for (model.name in names(model.varImp)) {
        imp.df <- model.varImp[[model.name]];
        var.importance[[model.name]] <- NA;
        var.importance[rownames(imp.df), model.name] <- imp.df$Overall;
    }
    if (include.ranks) {
        ranks <- as.data.frame(lapply(var.importance[, -1] * -1, rank));
        names(ranks) <- paste0(names(ranks), '.rank');
        ranks$mean.rank <- rowMeans(ranks);
        var.importance <- cbind.data.frame(var.importance, ranks);
    }
    var.importance;
}

#' Creates a summary data frame from the given (nested models)
#' Uses Youden's J statistic to find optimal threshold
#'
#' @param models the models with first level = target and second level list the different models for that target
#' @param models.roc optional ROC objects corresponding to the models
#'
#' @return
#' @export
summarize.models <- function(models, models.roc = NULL) {
    if (is.null(models.roc)) {
        models.roc <- lapply(models, function(x) {
            lapply(x, function(m) {
                bestPreds <- with(m, merge(pred, bestTune));
                pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
            })
        })
    }

    thresholds.info.list <- lapply(names(models), function(x) {
        target.list <- models[[x]]

        res <- lapply(names(target.list), function(y) {
            m <- target.list[[y]]
            m.roc <- models.roc[[x]][[y]]
            threshold.roc <- as.numeric(coords(m.roc, 'best', transpose = TRUE)[1])

            thres.res <- threshold.summary.stats(m, threshold.roc, per.fold = FALSE);
            thres.res$model <- y;
            thres.res$target <- x;
            thres.res$threshold <- threshold.roc;

            thres.res;
        })

        do.call(rbind.data.frame, res);
    })

    do.call(rbind.data.frame, thresholds.info.list);
}

#' Summarize the sequence models
#'
#' @param models
#' @param models.roc
#'
#' @return
#' @export
summarize.seq.models <- function(models, models.roc = NULL) {
    if (is.null(models.roc)) {
        models.roc <- lapply(models, function(m) {
            bestPreds <- with(m, merge(pred, bestTune));
            pROC::roc(predictor = bestPreds$yes, response = bestPreds$obs, direction = '<', levels = c('no', 'yes'));
        })
    }

    thresholds.info.list <- lapply(seq_along(models), function(x) {
        m <- models[[x]]
        m.roc <- models.roc[[x]]
        threshold.roc <- as.numeric(coords(m.roc, 'best', transpose = TRUE)[1])

        thres.res <- threshold.summary.stats(m, threshold.roc, per.fold = FALSE);
        thres.res$group <- x;
        thres.res$threshold <- threshold.roc;
        thres.res
    })

    do.call(rbind.data.frame, thresholds.info.list);
}
