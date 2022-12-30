#' Plot results of 'analysis_estimation'.
#'
#' @param results Results of 'analysis_estimation', including null distribution and bootstrap.
#'  Must be a list with at least the components: defs, table_null, boot_out.
#' @param sampsize Sample size of estimation data set.
#' @param alpha Size of confidence interval.
#' @param m_factor Fractional power for calculating resample size.
#' @param A1_length Number of levels of A1.
#' @param A2_length Number of levels of A2.
#' @param plot_labels Vector of (shortened) protected group names for labeling counterfactual error
#'  rate plots. Must have length(A1)*length(A2) components.
#' @param plot_values Vector of numeric ggplot2 shape scale values
#'  for counterfactual error rate plots. Must have length(A1)*length(A2) components.
#' @param plot_colors Vector of colors (string name or hex) for counterfactual error
#'  rate plots. Must have length(A1)*length(A2) components.
#'
#' @returns List of plots with the following components:
#'  * cfpr: Counterfactual false positive rates by group.
#'  * cfnr: Counterfactual false negative rates by group.
#'  * metrics_pos: Positive unfairness metrics.
#'  * metrics_neg: Negative unfairness metrics.
#'  * null_dist: Null distribution plots with estimated metrics and confidence intervals.
#'
#' @export

get_plots <- function(results, sampsize, alpha, m_factor, A1_length, A2_length,
                      plot_labels, plot_values, plot_colors) {
  ##################
  # Prepare results
  ##################
  # Named unfairness metrics
  est_named <- names_defs(results$defs, A1_length = A1_length, A2_length = A2_length)

  # Rescaled bootstrap results
  bs_rescaled <- get_bs_rescaled(bs_table = results$boot_out$ipw, defs = results$defs, sampsize = sampsize,
                                 A1_length = A1_length, A2_length = A2_length, m_factor = m_factor)

  # Confidence intervals
  ## t-interval
  ci_t_avg_neg <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_avg_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_avg_pos <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_avg_pos_new", sampsize = sampsize, alpha = alpha, m_factor = m_factor)

  ci_t_max_neg <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_max_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_max_pos <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_max_pos_new", sampsize = sampsize, alpha = alpha, m_factor = m_factor)

  ci_t_var_neg <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_var_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_var_pos <- ci_tint(bs_rescaled, est_named, parameter = "cdelta_var_pos_new", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ## t-interval, truncated
  ci_trunc_avg_neg <- ci_trunc(ci_t_avg_neg, type = "tint")
  ci_trunc_avg_pos <- ci_trunc(ci_t_avg_pos, type = "tint")

  ci_trunc_max_neg <- ci_trunc(ci_t_max_neg, type = "tint")
  ci_trunc_max_pos <- ci_trunc(ci_t_max_pos, type = "tint")

  ci_trunc_var_neg <- ci_trunc(ci_t_var_neg, type = "tint")
  ci_trunc_var_pos <- ci_trunc(ci_t_var_pos, type = "tint")

  ci_trunc_df <- list("cdelta_avg_neg" = ci_trunc_avg_neg, "cdelta_avg_pos_new" = ci_trunc_avg_pos,
                      "cdelta_max_neg" = ci_trunc_max_neg, "cdelta_max_pos_new" = ci_trunc_max_pos,
                      "cdelta_var_neg" = ci_trunc_var_neg, "cdelta_var_pos_new" = ci_trunc_var_pos) %>% dplyr::bind_rows(.id = "stat")
  rownames(ci_trunc_df) <- NULL
  ### For counterfactual error rates
  results_cfpr <- as.list(dplyr::select(est_named, (dplyr::contains("cfpr") & !dplyr::contains("marg"))))
  results_cfnr <- as.list(dplyr::select(est_named, (dplyr::contains("cfnr") & !dplyr::contains("marg"))))

  ci_trunc_cfpr <- lapply(names(results_cfpr), function(n) {
    ci_trunc(ci_tint(bs_rescaled, est_named, parameter = n, sampsize = sampsize, alpha = alpha, m_factor = m_factor), type = "tint")
  })
  names(ci_trunc_cfpr) <- names(results_cfpr)
  ci_trunc_cfpr <- ci_trunc_cfpr %>% dplyr::bind_rows(.id = "stat")

  ci_trunc_cfnr <- lapply(names(results_cfnr), function(n) {
    ci_trunc(ci_tint(bs_rescaled, est_named, parameter = n, sampsize = sampsize, alpha = alpha, m_factor = m_factor), type = "tint")
  })
  names(ci_trunc_cfnr) <- names(results_cfnr)
  ci_trunc_cfnr <- ci_trunc_cfnr %>% dplyr::bind_rows(.id = "stat")

  # Named null distribution results
  table_null <- names_df(results$table_null, A1_length = A1_length, A2_length = A2_length)
  # U-value table
  table_uval <- data.frame(
    cdelta_avg_neg = mean(table_null$cdelta_avg_neg < est_named$cdelta_avg_neg, na.rm = T),
    cdelta_avg_pos_new = mean(table_null$cdelta_avg_pos_new < est_named$cdelta_avg_pos_new, na.rm = T),
    cdelta_max_neg = mean(table_null$cdelta_max_neg < est_named$cdelta_max_neg, na.rm = T),
    cdelta_max_pos_new = mean(table_null$cdelta_max_pos_new < est_named$cdelta_max_pos_new, na.rm = T),
    cdelta_var_neg = mean(table_null$cdelta_var_neg < est_named$cdelta_var_neg, na.rm = T),
    cdelta_var_pos_new = mean(table_null$cdelta_var_pos_new < est_named$cdelta_var_pos_new, na.rm = T)
  ) %>% cond_round_3()
  ## Reversed (1-) u-value table
  table_uval_rev <- table_uval %>% dplyr::mutate(dplyr::across(dplyr::everything(), function(x) {1-as.numeric(x)}))

  # Estimation table with confidence intervals
  ## Put all CI information in same table
  allcis_trunc <- dplyr::bind_rows(list(ci_trunc_df, ci_trunc_cfpr, ci_trunc_cfnr))
  rownames(allcis_trunc) <- NULL
  ## Add CIs to estimation table
  est_summaries <- est_named %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "stat") %>%
    dplyr::mutate(sign = dplyr::case_when(
      stringr::str_detect(.data$stat, "cfpr_marg") ~ "cfpr_marg",
      stringr::str_detect(.data$stat, "cfnr_marg") ~ "cfnr_marg",
      stringr::str_detect(.data$stat, "cfpr") ~ "cfpr",
      stringr::str_detect(.data$stat, "cfnr") ~ "cfnr",
      stringr::str_detect(.data$stat, "fpr") ~ "fpr",
      stringr::str_detect(.data$stat, "fnr") ~ "fnr",
      .data$stat %in% c("cdelta_avg_neg", "cdelta_max_neg", "cdelta_var_neg") ~ "aggregate_neg",
      .data$stat %in% c("cdelta_avg_pos_new", "cdelta_max_pos_new", "cdelta_var_pos_new") ~ "aggregate_pos",
      TRUE ~ "other"
    )) %>% dplyr::left_join(allcis_trunc, by = "stat")

  #########
  # Plots
  #########
  # Counterfactual error rate plots
  p_cfpr <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfpr"), ggplot2::aes(x = "cFPR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFPR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_color_discrete(type = plot_colors, labels = plot_labels) +
    ggplot2::scale_shape_manual(values = plot_values, labels = plot_labels)

  p_cfnr <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfnr"), ggplot2::aes(x = "cFNR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFNR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_color_discrete(type = plot_colors, labels = plot_labels) +
    ggplot2::scale_shape_manual(values = plot_values, labels = plot_labels)

  # Unfairness metric plots
  p_metrics_neg <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "aggregate_neg"), ggplot2::aes(x = "Unfairness metric", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Negative metrics (Estimate, 95% CI)", title = NULL, shape = "Metric", color = "Metric") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_color_discrete(
      type = c("black", "red", "dodgerblue"), labels = c("Average", "Maximum", "Variational")) +
    ggplot2::scale_shape_manual(
      values = c(15, 16, 17), labels = c("Average", "Maximum", "Variational"))

  p_metrics_pos <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "aggregate_pos"), ggplot2::aes(x = "Unfairness metric", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Positive metrics (Estimate, 95% CI)", title = NULL, shape = "Metric", color = "Metric") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_color_discrete(
      type = c("black", "red", "dodgerblue"), labels = c("Average", "Maximum", "Variational")) +
    ggplot2::scale_shape_manual(
      values = c(15, 16, 17), labels = c("Average", "Maximum", "Variational"))

  # Null distribution/estimation plots
  ## Average
  p_null_neg_avg <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_avg_neg)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_avg_neg"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_avg_neg, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Average (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_neg)),
               x = 0.55, y = 1, hjust = "right", size = 5)

  p_null_pos_avg <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_avg_pos_new)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_avg_pos_new"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_avg_pos, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Average (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_pos_new)),
               x = 0.17, y = 5, hjust = "right", size = 5)
  ## Maximum
  p_null_neg_max <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_max_neg)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_max_neg"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_max_neg, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Maximum (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_neg)),
               x = 1, y = .5, hjust = "right", size = 5)

  p_null_pos_max <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_max_pos_new)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_max_pos_new"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_max_pos, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Maximum (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_pos_new)),
               x = 0.3, y = 2.5, hjust = "right", size = 5)
  ## Variational
  p_null_neg_var <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_var_neg)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_var_neg"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_var_neg, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Variational (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", round(.data$cdelta_var_neg, 3))),
               x = 0.085, y = 10, hjust = "right", size = 5)

  p_null_pos_var <- ggplot2::ggplot(table_null, ggplot2::aes(x = .data$cdelta_var_pos_new)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(data = dplyr::filter(est_summaries, .data$stat == "cdelta_var_pos_new"), ggplot2::aes(x = .data$value, y = 0), size = 3) +
    ggplot2::geom_segment(data = ci_trunc_var_pos, ggplot2::aes(x = .data$low_trans, xend = .data$high_trans), y = 0, yend = 0, size = 1) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Variational (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12)) +
    ggplot2::geom_label(data = table_uval_rev,
               mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_var_pos_new)),
               x = 0.011, y = 150, hjust = "right", size = 5)
  ## All in single plot
  p_null_all <- egg::ggarrange(plots = list(p_null_neg_avg, p_null_pos_avg, p_null_neg_max,
                                            p_null_pos_max, p_null_neg_var, p_null_pos_var),
                               nrow = 3, draw = F)

  # Return plots
  return(list(cfpr = p_cfpr, cfnr = p_cfnr, metrics_pos = p_metrics_pos,
              metrics_neg = p_metrics_neg, null_dist = p_null_all))
}
