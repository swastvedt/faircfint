% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_est}
\alias{data_est}
\title{Simulated data for demonstration of intersectional, counterfactual unfairness metrics.}
\format{
\subsection{\code{data_est}}{

A data frame with 5,000 rows and 14 columns:
\describe{
\item{A1, A2}{Protected characteristics (binary)}
\item{X.1, X.2, X.3, X.4}{Covariates (numeric)}
\item{Y0, Y1}{Potential outcomes (binary)}
\item{Y}{Observed outcome (binary)}
\item{D}{Treatment (binary)}
\item{S}{Binary risk prediction (using cut-off of 0.5)}
\item{S_prob}{Continuous risk prediction}
...
}
}
}
\source{
Generated via simulation.
}
\usage{
data_est
}
\description{
Simulated data for demonstration of intersectional, counterfactual unfairness metrics.
}
\keyword{datasets}
