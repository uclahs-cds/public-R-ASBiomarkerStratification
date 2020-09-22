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
#'
#' @examples
#' camel.to.spaces(c("BiopsyUpgraded", "ProgressedToTreatment"))
camel.to.spaces <- function(x) {
    matches <- gsub("(?!^)([[:upper:]])", " \\1", x, perl = TRUE)
    split.string <- strsplit(matches, " ")
    unlist(lapply(split.string, paste0, collapse = " "))
}
