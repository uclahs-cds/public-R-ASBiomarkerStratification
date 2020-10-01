#' Title
#'
#' @param biodb
#' @param y
#' @param x
#' @param cond
#'
#' @return
#' @export
#'
#' @examples
demographics.boxplot <- function(biodb, y, x, cond = NULL, ...) {
    box.formula <- paste(y, x, sep = " ~ ")
    if(!is.null(cond)) {
        box.formula <- paste(box.formula, cond, sep = " | ");
        }
    ylab.label <- attr(biodb[, y], 'label');
    if(is.null(ylab.label)) {
        ylab.label <- y;
        }
    biodb[, cond] <- as.factor(biodb[, cond]);

    create.boxplot(
        formula = as.formula(box.formula),
        data = biodb,
        add.stripplot = TRUE,
        xlab.label = camel.to.spaces(x),
        ylab.label = ylab.label,
        # ylimits = c(0, max(biodb[, y] + 10)),
        ...
    );
}

#' Title
#'
#' @param models.roc
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
roc.pr.plot <- function(models.roc, ...) {
    opar <- par(pty="s", mfrow=c(1,2))
    on.exit(par(opar))

    plot(
        x = models.roc[[1]],
        las = 1,
        col = default.colours(2)[1],
        main = 'ROC Curve',
        ylim = c(0, 1),
        # pROC resets mar and mpg..
        # https://github.com/xrobin/pROC/blob/master/R/plot.roc.R#L150
        mar = par('mar'),
        mgp = par('mgp'),
        print.thres = TRUE,
        # print.thres.adj=c(1, -.6),
        print.thres.cex=.75,
        print.thres.pattern = "%.2f",
        ...
    )
    if(2 == length(models.roc)) {
        plot(models.roc[[2]], las = 1, add = TRUE, col = default.colours(2)[2], xlim = c(1, 0), print.thres = TRUE,
             print.thres.cex=.75, print.thres.pattern = "%.2f", ...);
    }

    legend("bottomright", names(models.roc), cex = 0.8, col = default.colours(2), lwd = c(2, 2), inset = 0.02);

    # Plot PR Curves
    # Plot the Precision-Recall curve
    plot(precision ~ recall,
         coords(models.roc[[1]], "all", ret = c("recall", "precision"), transpose = FALSE),
         type="l", las = 1, col = default.colours(2)[1], lwd = 2, main = 'Precision-Recall Curve', ylim = c(0, 1), ...);

    if(2 == length(models.roc)) {
        lines(precision ~ recall,
              coords(models.roc[[2]], "all", ret = c("recall", "precision"), transpose = FALSE),
              type="l", las = 1, col = default.colours(2)[2], lwd = 2, ...);
        }
}
