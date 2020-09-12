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
