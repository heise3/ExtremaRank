\name{extremarank_methods}
\alias{summary.extremarank_result}
\alias{print.summary.extremarank_result}
\alias{as.data.frame.extremarank_result}
\alias{plot.extremarank_result}
\alias{plot.extremarank_curve}
\title{Summarize and Plot Deletion Audits}
\description{Summaries expose coverage and candidate status. Plots distinguish exact native certificates from observed model-refit results.}
\usage{
\method{summary}{extremarank_result}(object, ...)
\method{print}{summary.extremarank_result}(x, ...)
\method{as.data.frame}{extremarank_result}(x, row.names = NULL, optional = FALSE, ...)
\method{plot}{extremarank_result}(x, type = c("status", "influence"), ...)
\method{plot}{extremarank_curve}(x, ...)
}
\arguments{
\item{object, x}{An extremarank_result or, for the curve plot, an extremarank_curve.}
\item{type}{Candidate membership counts or the Top-K Jaccard change of actually evaluated deletions. Missing/failed scenarios are not assigned zero influence.}
\item{row.names, optional}{Data-frame generic arguments; result feature rows are preserved.}
\item{...}{Additional plotting arguments, or unused generic arguments.}
}
\value{summary returns a summary.extremarank_result with coverage, timings and status counts. as.data.frame returns candidate results or the fit-status table. Plot and print methods return the input invisibly.}
\seealso{\code{\link{extremarank}}, \code{\link{robustness_curve}}}
