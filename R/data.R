#' Simulated data for demonstration of intersectional, counterfactual unfairness metrics.
#'
#' @format ## `ex_data_estimation`
#' A data frame with 5,000 rows and 14 columns:
#' \describe{
#'   \item{A1, A2}{Protected characteristics (binary)}
#'   \item{X.1, X.2, X.3, X.4}{Covariates (numeric)}
#'   \item{Y}{Observed outcome (binary)}
#'   \item{D}{Treatment (binary)}
#'   \item{S}{Binary risk prediction (using cut-off of 0.5)}
#'   \item{S_prob}{Continuous risk prediction}
#'   ...
#' }
#' @source Generated via simulation.
"ex_data_estimation"

#' Simulated data for demonstration of small subgroup unfairness metrics.
#'
#' @format ## `ex_data_small`
#' A data frame with 5,000 rows and 25 columns:
#' \describe{
#'   \item{A1, A2}{Protected characteristics}
#'   \item{X.1, X.2, X.3, X.4}{Covariates for propensity score model (numeric)}
#'   \item{X_pa.1, ..., X_pa.7}{Covariates for P(A=a) model (numeric)}
#'   \item{X_outcome.1, ..., X_outcome.8}{Covariates for outcome models (numeric)}
#'   \item{Y}{Observed outcome (binary)}
#'   \item{D}{Treatment (binary)}
#'   \item{S}{Binary risk prediction (using cut-off of 0.5)}
#'   \item{S_prob}{Continuous risk prediction}
#'   ...
#' }
#' @source Generated via simulation.
"ex_data_small"

#' Simulated 'external data' for demonstration of small subgroup unfairness metrics with data borrowing.
#'
#' @format ## `ex_data_external`
#' A data frame with 5,000 rows and 25 columns:
#' \describe{
#'   \item{A1, A2}{Protected characteristics}
#'   \item{X_pa.1, ..., X_pa.5}{Covariates for external P(A=a) model (numeric)}
#'   ...
#' }
#' @source Generated via simulation.
"ex_data_external"
