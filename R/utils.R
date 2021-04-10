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
#' camel.to.spaces(c('BiopsyUpgraded', 'ProgressedToTreatment'))
camel.to.spaces <- function(x, replace = ' ') {
    matches <- gsub('(?!^)([[:upper:]])', ' \\1', x, perl = TRUE);
    split.string <- strsplit(matches, ' ');
    unlist(lapply(split.string, paste0, collapse = replace));
  }

# Moving average code from forecast
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
        stop('order must be an integer')
        };
    if (order %% 2 == 0 && centre) {
        w <- c(0.5, rep(1, order - 1), 0.5) / order
        };
    else {
        w <- rep(1, order) / order
        };
    return(stats::filter(x, w));
    }

label <- function(x) {
    attr(x, 'label');
    }

#' Get the label or the name of a variable or data frame. If x is a data frame, then return the label or column name.
#'
#' @param x the variable or data frame we want to get either the label or name from.
#'
#' @export
#' @examples
#' x <- data.frame(x = 1:4, y = 5:8)
#' attr(x$x, 'label') <- 'Label 1'
#' label.or.name(x$x)
#' #> "Label 1"
#' label.or.name(x$y)
#' #> x$y
#' label.or.name(x)
#' #> "Label 1" "y"
label.or.name <- function(x) {
    mc <- match.call()
    if (typeof(x) == 'list') {
        res <- names(x)
        # Nulls are dropped
        labs <- lapply(x, label)
        cond <- ! unlist(lapply(labs, is.null))
        res[cond] <- unlist(labs)
        res
        }
    else {
        lab <- label(x)
        if (is.null(lab)) mc$x
        else lab
        }
    }
