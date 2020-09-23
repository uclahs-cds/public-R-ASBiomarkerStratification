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

