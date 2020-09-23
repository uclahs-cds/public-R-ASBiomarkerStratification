#' Render a simple binary table1 table to
#'
#' @param tb1
#' @param group
#'
#' @return
#' @export
#'
#' @examples
table1.to.kable <- function(tb1, group) {
    rt <- rvest::html_table(xml2::read_html(tb1))[[1]]
    group.index <- which(rt[,2] == "");
    group.names <- rt[group.index, 1];

    # Remove the group rows
    rt.nogroup <- rt[-group.index, ];

    k1 <- kable(rt.nogroup, align = c('l', 'c', 'c', 'c'), row.names = FALSE);

    group.index.aligned <- group.index - seq(0,length(group.index) - 1)
    group.index.end <- c(group.index.aligned, nrow(rt.nogroup))
    for(i in 1:length(group.names)) {
        k1 <- pack_rows(k1, group.names[i], group.index.end[i], group.index.end[i + 1], indent = FALSE)
    }

    add_indent(k1, 1:nrow(rt.nogroup))
}
