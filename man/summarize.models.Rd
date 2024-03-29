% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval.models.R
\name{summarize.models}
\alias{summarize.models}
\title{Creates a summary data frame from the given (nested models)
Uses Youden's J statistic to find optimal threshold}
\usage{
summarize.models(models, models.roc = NULL)
}
\arguments{
\item{models}{the models with first level = target and second level list the different models for that target}

\item{models.roc}{optional ROC objects corresponding to the models}
}
\value{

}
\description{
Creates a summary data frame from the given (nested models)
Uses Youden's J statistic to find optimal threshold
}
