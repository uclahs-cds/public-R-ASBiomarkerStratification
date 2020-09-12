#' Converts event times and end of study times into survival::Surv format
#'
#' @param event.times The times of the events. Expects NA if event did not occur.
#' @param end.times The last follow up time.
#'
#' @return a list with keys `time` and `event`
#' @export
#'
#' @examples
#' x <- c(NA, 1, NA, 2)
#' y <- c(3, 4, 4.5, 7)
#' do.call(survival::Surv, surv.format(x, y))
surv.format <- function(event.times, end.times) {
    time <- event.times
    censored <- is.na(event.times)
    time[censored] <- end.times[censored]
    event <- ifelse(!censored, 1, 0)
    # If we are missing both the event time and end time then event is NA
    event[is.na(end.times) & censored] <- NA
    list(
        time = time,
        event = event
        )
}
