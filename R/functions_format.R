#' Extract names of variables with non-zero coefficients after cv.glmnet.
#'
#' @param model cv.glmnet model object.
#' @param lambda Optional lambda value for selecting coefficients. Defaults to lambda.1se.
#' @param xvars Variable names for all predictors in cv.glmnet object.
#' @param outcome Outcome variable name (for formula)
#'
#' @returns A list with the following components:
#'  * coef: Dataframe with non-zero coefficients
#'  * vars_names: Names of variables with non-zero coefficients
#'  * formula: Formula with selected variables

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

#' Get vector of names for output of 'get_defs_analysis' function.
#'
#' @param A1_length Length of 'A1' protected characteristic variable.
#' @param A2_length Length of 'A2' protected characteristic variable.
#'
#' @returns Vector of names.

defs_names <- function(A1_length, A2_length) {
  A1A2_length <- A1_length*A2_length
  c(paste("fpr", as.character(seq(1:A1A2_length)), sep = "_"),
    paste("fnr", as.character(seq(1:A1A2_length)), sep = "_"),
    paste("cfpr", as.character(seq(1:A1A2_length)), sep = "_"),
    paste("cfnr", as.character(seq(1:A1A2_length)), sep = "_"),
    paste("cfpr_marg_A1", as.character(seq(1:A1_length)), sep = "_"),
    paste("cfpr_marg_A2", as.character(seq(1:A2_length)), sep = "_"),
    paste("cfnr_marg_A1", as.character(seq(1:A1_length)), sep = "_"),
    paste("cfnr_marg_A2", as.character(seq(1:A2_length)), sep = "_"),
    "cdelta_avg_neg", "cdelta_max_neg", "cdelta_quant_neg", "cdelta_var_neg", "odelta_int_neg", "cdelta_marg_neg",
    "cdelta_avg_pos_new", "cdelta_max_pos_new", "cdelta_quant_pos_new", "cdelta_var_pos_new", "odelta_int_pos", "cdelta_marg_pos_new"
  )
}

#' Assign column names to output of 'get_defs_analysis' or other list result object.
#'
#' @inheritParams defs_names
#' @param defs List of definitions, output by 'get_defs_analysis' function.
#'
#' @returns Dataframe of unfairness metrics with named columns.

names_defs <- function(defs, A1_length, A2_length) {
  df <- as.data.frame(t(defs))
  colnames(df) <- defs_names(A1_length, A2_length)
  return(df)
}

#' Assign column names to bootstrap results or other dataframe result object.
#'
#' @inheritParams defs_names
#' @param bs_table Dataframe of bootstrap results.
#'
#' @returns Dataframe of bootstrap results with named columns.

names_df <- function(bs_table, A1_length, A2_length) {
  df <- as.data.frame(bs_table)
  colnames(df) <- defs_names(A1_length, A2_length)
  return(df)
}

#' Create rescaled bootstrap estimate table from results of bs_rescaled_analysis.
#'
#' @inheritParams names_df
#' @inheritParams names_defs
#' @param sampsize Sample size of estimation data set.
#' @param m_factor Fractional power for calculating resample size (default 0.75).
#'
#' @returns Dataframe of rescaled bootstrap estimates.

get_bs_rescaled <- function(bs_table, defs, sampsize, A1_length, A2_length,
                            m_factor) {
  bs_table <- names_df(bs_table, A1_length = A1_length, A2_length = A2_length)
  est_vals <- names_defs(defs, A1_length = A1_length, A2_length = A2_length)[1,]
  m <- sampsize^m_factor
  sqrt(m)*sweep(bs_table, MARGIN = 2, STATS = unlist(est_vals), FUN = "-")
}

#' Calculate normal approximation confidence interval.
#'
#' @param bs_rescaled Rescaled bootstrap estimate table (result of get_bs_rescaled).
#' @param est_named Named 1-row dataframe of raw metric estimates (result of names_defs).
#' @param parameter Column name of the metric being estimated (string).
#' @param sampsize Sample size of estimation data set.
#' @param alpha Size of (1-alpha) confidence interval.
#'
#' @returns Dataframe with point estimate, SE, variance, 1-alpha confidence interval.

ci_norm <- function(bs_rescaled, est_named, parameter, sampsize, alpha) {
  data.frame(
    var_est = stats::var(bs_rescaled[[parameter]], na.rm = T),
    point_est = est_named[[parameter]]
  ) %>%
  dplyr::mutate(
    se_est = sqrt(.data$var_est/sampsize),
    ci_low = .data$point_est - stats::qnorm(1-alpha/2)*.data$se_est,
    ci_high = .data$point_est + stats::qnorm(1-alpha/2)*.data$se_est
  )
}

#' Calculate t-interval.
#'
#' @inheritParams ci_norm
#' @param m_factor Fractional power for calculating resample size (default 0.75).
#'
#' @returns Dataframe with point estimate, SE, variance, 1-alpha confidence interval.

ci_tint <- function(bs_rescaled, est_named, parameter, sampsize, alpha, m_factor) {
  se_table <- ci_norm(bs_rescaled, est_named, parameter, sampsize, alpha)
  data_temp <- bs_rescaled %>%
    dplyr::mutate(t_value = bs_rescaled[[parameter]]/(sqrt(sampsize^m_factor)*se_table$se_est))
  data.frame(
    point_est = se_table$point_est,
    se_est = se_table$se_est,
    low_trans = se_table$point_est - se_table$se_est*stats::quantile(data_temp$t_value, 1-alpha/2, na.rm = T),
    high_trans = se_table$point_est - se_table$se_est*stats::quantile(data_temp$t_value, alpha/2, na.rm = T)
  )
}

#' Truncate confidence interval.
#'
#' @param ci_result Confidence interval dataframe (result of ci_norm or ci_tint).
#' @param type Type of confidence interval ('norm' or 'tint')
#'
#' @returns Confidence interval dataframe with truncated interval.

ci_trunc <- function(ci_result, type) {
  if(type == "norm") {
    ci_result %>%
      dplyr::mutate(ci_low = max(0,.data$ci_low), ci_high = min(1,.data$ci_high))
  } else if(type == "tint") {
    ci_result %>%
      dplyr::mutate(low_trans = max(0,.data$low_trans), high_trans = min(1,.data$high_trans))
  }
}
