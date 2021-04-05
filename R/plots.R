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
    box.formula <- paste(y, x, sep = ' ~ ')
    if (!is.null(cond)) {
        box.formula <- paste(box.formula, cond, sep = ' | ');
        }
    ylab.label <- attr(biodb[, y], 'label');
    if (is.null(ylab.label)) {
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
    opar <- par(pty = 's', mfrow = c(1,2), mar = 0.1 + c(8, 6, 4, 2), cex = 1.3, cex.axis = 1.1, cex.lab = 1.5)
    on.exit(par(opar))

    models.n <- length(models.roc);

    plot(
        x = models.roc[[1]],
        las = 1,
        col = BoutrosLab.plotting.general::default.colours(models.n)[1],
        main = 'ROC Curve',
        ylim = c(0, 1),
        # pROC resets mar and mpg..
        # https://github.com/xrobin/pROC/blob/master/R/plot.roc.R#L150
        mar = par('mar'),
        mgp = par('mgp'),
        print.thres = TRUE,
        # print.thres.adj=c(1, -.6),
        #print.thres.cex=.75,
        print.thres.pattern = '',
        print.thres.col = BoutrosLab.plotting.general::default.colours(models.n)[1],
        cex.axis = par('cex.axis'),
        cex.lab = par('cex.lab'),
        #cex = par('cex'),
        # print.auc=TRUE,
        ...
    )
    for (i in 2:models.n) {
        plot(
            models.roc[[i]],
            las = 1,
            add = TRUE,
            col = BoutrosLab.plotting.general::default.colours(models.n)[i],
            xlim = c(1, 0),
            print.thres = TRUE,
            # print.thres.cex=.75,
            print.thres.pattern = '', # '%.2f',
            print.thres.col = BoutrosLab.plotting.general::default.colours(models.n)[i],
            cex.axis = par('cex.axis'),
            cex.lab = par('cex.lab'),
            # cex = par('cex'),
            # print.auc=TRUE,
            # print.auc.adj = c(0, i),
            ...
        );
    }

    # Plot PR Curves
    # Plot the Precision-Recall curve
    plot(precision ~ recall,
         coords(models.roc[[1]], 'all', ret = c('recall', 'precision'), transpose = FALSE),
         type = 'l', las = 1, col = BoutrosLab.plotting.general::default.colours(models.n)[1], lwd = 2, main = 'Precision-Recall Curve', ylim = c(0, 1),
         ylab = 'Precision',
         xlab = 'Recall',
         cex.axis = par('cex.axis'),
         cex.lab = par('cex.lab'),
         cex = par('cex'), ...);

    for (i in 2:models.n) {
        lines(precision ~ recall,
              coords(models.roc[[i]], 'all', ret = c('recall', 'precision'), transpose = FALSE),
              type = 'l', las = 1, col = BoutrosLab.plotting.general::default.colours(models.n)[i], lwd = 2,
              ylab = 'Precision',
              xlab = 'Recall',
              cex.axis = par('cex.axis'),
              cex.lab = par('cex.lab'),
              cex = par('cex'), ...);
    }

    # models.auc <- unlist(lapply(models.roc, `[[`, i = 'auc'));

    auc.ci.text <- lapply(models.roc, function(m) {
        auc.ci <- ci(m$auc);
        sprintf('%.2f (%.2f, %.2f) 95%% CI', auc.ci[2], auc.ci[1], auc.ci[3])
        });

    legend.text <- paste0(names(models.roc), ': AUC = ', auc.ci.text);

    legend(-1.8, -0.5, legend.text, cex = 1.1, col = BoutrosLab.plotting.general::default.colours(models.n), lwd = c(2, 2), inset = 0.02, xpd = 'NA');
}

roc.pr.multiple.plot <- function(models.roc, ...) {
    opar <- par(pty = 's', mfrow = c(1,2), mar = 0.1 + c(8, 6, 4, 2), cex = 1.3, cex.axis = 1.1, cex.lab = 1.5)
    on.exit(par(opar))

    models.n <- length(models.roc);

    plot(
        x = models.roc[[1]],
        las = 1,
        col = BoutrosLab.plotting.general::default.colours(models.n)[1],
        main = 'ROC Curve',
        ylim = c(0, 1),
        # pROC resets mar and mpg..
        # https://github.com/xrobin/pROC/blob/master/R/plot.roc.R#L150
        mar = par('mar'),
        mgp = par('mgp'),
        print.thres = TRUE,
        # print.thres.adj=c(1, -.6),
        #print.thres.cex=.75,
        print.thres.pattern = '',
        print.thres.col = BoutrosLab.plotting.general::default.colours(models.n)[1],
        cex.axis = par('cex.axis'),
        cex.lab = par('cex.lab'),
        #cex = par('cex'),
        # print.auc=TRUE,
        ...
    )
    for (i in 2:models.n) {
        plot(
            models.roc[[i]],
            las = 1,
            add = TRUE,
            col = BoutrosLab.plotting.general::default.colours(models.n)[i],
            xlim = c(1, 0),
            print.thres = TRUE,
            # print.thres.cex=.75,
            print.thres.pattern = '', # '%.2f',
            print.thres.col = BoutrosLab.plotting.general::default.colours(models.n)[i],
            cex.axis = par('cex.axis'),
            cex.lab = par('cex.lab'),
            # cex = par('cex'),
            # print.auc=TRUE,
            # print.auc.adj = c(0, i),
            ...
        );
    }

    # Plot PR Curves
    # Plot the Precision-Recall curve
    plot(precision ~ recall,
         coords(models.roc[[1]], 'all', ret = c('recall', 'precision'), transpose = FALSE),
         type = 'l', las = 1, col = BoutrosLab.plotting.general::default.colours(models.n)[1], lwd = 2, main = 'Precision-Recall Curve', ylim = c(0, 1),
         ylab = 'Precision',
         xlab = 'Recall',
         cex.axis = par('cex.axis'),
         cex.lab = par('cex.lab'),
         cex = par('cex'), ...);

    for (i in 2:models.n) {
        lines(precision ~ recall,
              coords(models.roc[[i]], 'all', ret = c('recall', 'precision'), transpose = FALSE),
              type = 'l', las = 1, col = BoutrosLab.plotting.general::default.colours(models.n)[i], lwd = 2,
              ylab = 'Precision',
              xlab = 'Recall',
              cex.axis = par('cex.axis'),
              cex.lab = par('cex.lab'),
              cex = par('cex'), ...);
    }

    # models.auc <- unlist(lapply(models.roc, `[[`, i = 'auc'));

    auc.ci.text <- lapply(models.roc, function(m) {
        auc.ci <- ci(m$auc);
        sprintf('%.2f (%.2f, %.2f) 95%% CI', auc.ci[2], auc.ci[1], auc.ci[3])
    });

    legend.text <- paste0(names(models.roc), ': AUC = ', auc.ci.text);

    legend(-1.8, -0.5, legend.text, cex = 1.1, col = BoutrosLab.plotting.general::default.colours(models.n), lwd = c(2, 2), inset = 0.02, xpd = 'NA');
}
