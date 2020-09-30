preds.cost <- function(x, lev = NULL, cost.matrix = matrix(c(0,1,1,0), byrow = TRUE, nrow = 2), threshold) {
    if(is.null(lev)) {
        lev <- levels(x$obs)
    }
    preds <- x[[lev[2]]]
    truth <- x$obs
    class_preds <- as.factor(ifelse(preds > threshold, lev[2], lev[1]))
    levels(class_preds) <- lev
    conf.mat <- table(class_preds, truth)
    sum(conf.mat * cost.matrix)
}

#' Title
#'
#' @param model
#' @param cost.matrix
#' @param thresholds
#'
#' @return
#' @export
#'
#' @examples
cost.threshold.train <- function(model,
                              cost.matrix = matrix(c(0,1,1,0), byrow = TRUE, nrow = 2),
                              thresholds = 0.5) {
    stopifnot(class(model) == 'train');

    pred_best <- with(model,
                      merge(pred, bestTune))

    res <- lapply(thresholds, function(threshold) {
        grp_resample <- split(pred_best, pred_best$Resample)
        # Maybe add summary here?
        mean(unlist(lapply(grp_resample, preds.cost, threshold = threshold, cost.matrix = cost.matrix)))
    })
    # names(res) <- thresholds
    # If we are using summary
    # t(do.call(cbind, res))
    unlist(res)
}

#' Title
#'
#' @param model
#' @param cost.matrix
#' @param thresholds
#' @param smooth
#'
#' @return
#' @export
#'
#' @examples
optimal.threshold.train <- function(model,
                                    cost.matrix = matrix(c(0,1,1,0), byrow = TRUE, nrow = 2),
                                    thresholds = c(0.25, 0.5, 0.75),
                                    smooth = TRUE) {
    threshold.cost <- cost.threshold.train(model, thresholds = thresholds, cost.matrix = cost.matrix)
    if(smooth) {
        threshold.cost <- as.vector(moving.avg(threshold.cost, order = max(1, round(thresholds / 20))))
    }
    thresholds[which.min(threshold.cost)]
}
