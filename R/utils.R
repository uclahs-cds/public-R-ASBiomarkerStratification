#' Returns the group matches from a regular expression on a vector
#' @param x the vector we want to match on
#' @param pattern regular expression with groups
str.match <- function(x, pattern) {
    regmatches(x, regexec(pattern, x));
}


#' Convert camel case to string with spaces. Note does not work with multiple capitals
#' Idea from: https://stackoverflow.com/a/8407047/1351718
#'
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' camel.to.spaces(c("BiopsyUpgraded", "ProgressedToTreatment"))
camel.to.spaces <- function(x, replace = " ") {
    matches <- gsub("(?!^)([[:upper:]])", " \\1", x, perl = TRUE)
    split.string <- strsplit(matches, " ")
    unlist(lapply(split.string, paste0, collapse = replace))
}

# Moving average code from forecast
#' Title
#'
#' @param x
#' @param order
#' @param center
#'
#' @return
#' @export
#'
#' @examples
moving.avg <- function(x, order, center = TRUE) {
    if (abs(order - round(order)) > 1e-08) {
        stop("order must be an integer")
    }
    if (order%%2 == 0 && centre) {
        w <- c(0.5, rep(1, order - 1), 0.5)/order
    }
    else {
        w <- rep(1, order)/order
    }
    return(stats::filter(x, w))
}

# From fishmethods
combinevar <- function(xbar = NULL,s_squared = NULL,n = NULL) {
    if(length(xbar)!=length(s_squared)|length(xbar)!=length(n)|
       length(s_squared)!=length(n)) stop ("Vector lengths are different.")
    sum_of_squares <- sum((n-1)*s_squared + n*xbar^2)
    grand_mean <- sum(n*xbar)/sum(n)
    combined_var <- (sum_of_squares - sum(n)*grand_mean^2)/(sum(n)-1)
    return(c(grand_mean,combined_var))
}
