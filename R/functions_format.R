cond_round_3 <- function(table) {
  table %>% dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round, digits=3)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), function(x) {replace(x, x==0, "<0.001")}))
}
# utils::globalVariables("where")

# Extract names of variables with non-zero coefficients after cv.glmnet.
#
# @param model cv.glmnet model object.
# @param lambda Optional lambda value for selecting coefficients. Defaults to lambda.1se.
# @param xvars Variable names for all predictors in cv.glmnet object.
# @param outcome Outcome variable name (for formula)
#
# @returns A list with the following components:
#  * coef: Dataframe with non-zero coefficients
#  * vars_names: Names of variables with non-zero coefficients
#  * formula: Formula with selected variables

select_coef <- function(model, lambda = NULL, xvars, outcome) {
  if(!is.null(lambda)) {
    l = lambda
  } else {
    l = model$lambda.1se
  }
  coef <- data.frame(coef.name = dimnames(coef(model, s=l))[[1]],
                     coef.value = matrix(coef(model, s=l))) %>%
    dplyr::filter(.data$coef.value != 0) %>%
    dplyr::select(.data$coef.name)

  # Extract non-zero coefficients and paste into formula
  vars_regex <- paste(xvars, collapse="|")
  vars_names <- unique(stringr::str_extract_all(coef[,1], vars_regex, simplify = T))
  vars_names <- as.vector(vars_names[-1,])
  formula <- paste0(outcome, ' ~', paste(vars_names, collapse = " + " ))

  return(list(coef = coef, vars_names = vars_names, formula = formula))
}

# Create rescaled bootstrap estimate table from results of bs_rescaled_analysis.
#
# @param bs_table Dataframe of bootstrap results.
# @param est_vals List of definitions, output by 'get_defs_analysis' function.
# @param sampsize Sample size of estimation data set.
# @param m_factor Fractional power for calculating resample size.
#
# @returns Dataframe of rescaled bootstrap estimates.

get_bs_rescaled <- function(bs_table, est_vals, sampsize, m_factor) {
  m <- sampsize^m_factor
  sqrt(m)*sweep(bs_table, MARGIN = 2, STATS = unlist(est_vals), FUN = "-")
}

# Calculate normal approximation confidence interval.
#
# @param bs_rescaled Rescaled bootstrap estimate table (result of get_bs_rescaled).
# @param est_named Vector of metrics ('defs' element of'analysis_estimation' function output).
# @param parameter Name of the metric being estimated (string).
# @param sampsize Sample size of estimation data set.
# @param alpha Size of confidence interval.
#
# @returns Dataframe with point estimate, SE, variance, 1-alpha confidence interval.

ci_norm <- function(bs_rescaled, est_named, parameter, sampsize, alpha) {
  data.frame(
    var_est = stats::var(bs_rescaled[,parameter], na.rm = T),
    point_est = est_named[[parameter]]
  ) %>%
  dplyr::mutate(
    se_est = sqrt(.data$var_est/sampsize),
    ci_low = .data$point_est - stats::qnorm(1-alpha/2)*.data$se_est,
    ci_high = .data$point_est + stats::qnorm(1-alpha/2)*.data$se_est
  )
}

# Calculate t-interval.
#
# @inheritParams ci_norm
# @param m_factor Fractional power for calculating resample size.
#
# @returns Dataframe with point estimate, SE, variance, 1-alpha confidence interval.

ci_tint <- function(bs_rescaled, est_named, parameter, sampsize, alpha, m_factor) {
  se_table <- ci_norm(bs_rescaled, est_named, parameter, sampsize, alpha)
  data_temp <- as.data.frame(bs_rescaled) %>%
    dplyr::mutate(t_value = bs_rescaled[,parameter]/(sqrt(sampsize^m_factor)*se_table$se_est))
  data.frame(
    point_est = se_table$point_est,
    se_est = se_table$se_est,
    low_trans = se_table$point_est - se_table$se_est*stats::quantile(data_temp$t_value, 1-alpha/2, na.rm = T),
    high_trans = se_table$point_est - se_table$se_est*stats::quantile(data_temp$t_value, alpha/2, na.rm = T)
  )
}

# Truncate confidence interval.
#
# @param ci_result Confidence interval dataframe (result of ci_norm or ci_tint).
# @param type Type of confidence interval ('norm' or 'tint')
#
# @returns Confidence interval dataframe with truncated interval.

ci_trunc <- function(ci_result, type) {
  if(type == "norm") {
    ci_result %>%
      dplyr::mutate(ci_low = max(0,.data$ci_low), ci_high = min(1,.data$ci_high))
  } else if(type == "tint") {
    ci_result %>%
      dplyr::mutate(low_trans = max(0,.data$low_trans), high_trans = min(1,.data$high_trans))
  }
}
