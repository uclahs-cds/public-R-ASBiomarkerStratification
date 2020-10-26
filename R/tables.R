#' Render a simple binary table1 table to
#'
#' @param tb1
#' @param group
#' @param ... pass additional arguments into kable
#'
#' @return
#' @export
#'
#' @examples
table1.to.kable <- function(tb1, group, ...) {
    rt <- rvest::html_table(xml2::read_html(tb1))[[1]]
    group.index <- which(rt[,2] == "");
    group.names <- rt[group.index, 1];

    # Remove the group rows
    rt.nogroup <- rt[-group.index, ];

    k1 <- kable(rt.nogroup, align = c('l', 'c', 'c', 'c'), row.names = FALSE, booktabs = T, ...);

    group.index.aligned <- group.index - seq(0,length(group.index) - 1)
    group.index.end <- c(group.index.aligned, nrow(rt.nogroup))
    for(i in 1:length(group.names)) {
        k1 <- pack_rows(k1, group.names[i], group.index.end[i], group.index.end[i + 1], indent = FALSE)
    }

    add_indent(k1, 1:nrow(rt.nogroup))
}

euro.urology.table1 <- function(variables, data, cohort) {
    for(v in variables) {
        feature <- data[[v]]
        if(is.numeric(data[[v]])) {
            quartiles <- quantile(feature, probs = c(.25, .75, .5), na.rm = TRUE)
            print(sprintf("%.02f (%.02f, %.02f)", quartiles[1], quartiles[2], quartiles[3]))
        }
    }
}
