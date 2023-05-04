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
#' @param delta_uval Threshold for u-value test. Must be between 0 and 1.
#' @param uval_labels Vector of x and y coordinates for labels on u-value plots. Must have length
#'  12.
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
                      plot_labels, plot_values, plot_colors, delta_uval, uval_labels) {
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

  est_named_pivot <- est_named %>% dplyr::select(.data$cdelta_avg_neg, .data$cdelta_avg_pos_new,
                                                 .data$cdelta_max_neg, .data$cdelta_max_pos_new,
                                          .data$cdelta_var_neg, .data$cdelta_var_pos_new) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "stat", values_to = "value_obs")

  table_null_delta <- table_null %>% dplyr::select(.data$cdelta_avg_neg, .data$cdelta_avg_pos_new,
                                                   .data$cdelta_max_neg, .data$cdelta_max_pos_new,
                                            .data$cdelta_var_neg, .data$cdelta_var_pos_new) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "stat", values_to = "value_null") %>%
    dplyr::left_join(est_named_pivot, by = "stat") %>%
    dplyr::mutate(obs_minus_null = .data$value_obs - .data$value_null)

  # U-value table
  table_uval <- data.frame(
    cdelta_avg_neg = mean((est_named$cdelta_avg_neg - table_null$cdelta_avg_neg) > delta_uval, na.rm = T),
    cdelta_avg_pos_new = mean((est_named$cdelta_avg_pos_new - table_null$cdelta_avg_pos_new) > delta_uval, na.rm = T),
    cdelta_max_neg = mean((est_named$cdelta_max_neg - table_null$cdelta_max_neg) > delta_uval, na.rm = T),
    cdelta_max_pos_new = mean((est_named$cdelta_max_pos_new - table_null$cdelta_max_pos_new) > delta_uval, na.rm = T),
    cdelta_var_neg = mean((est_named$cdelta_var_neg - table_null$cdelta_var_neg) > delta_uval, na.rm = T),
    cdelta_var_pos_new = mean((est_named$cdelta_var_pos_new - table_null$cdelta_var_pos_new) > delta_uval, na.rm = T)
  ) %>% cond_round_3()

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
  p_null_neg_avg <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_avg_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Obs. - Null", y = "density", title = "Average (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_avgneg <- ggplot2::ggplot_build(p_null_neg_avg)$data[[1]]
  d_avgneg_sub <- subset(d_avgneg, .data$x > 0.1)
  if(nrow(d_avgneg_sub) > 0) {
    p_null_neg_avg <- p_null_neg_avg +
      ggplot2::geom_area(data = d_avgneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_neg)),
                          x = uval_labels[1], y = uval_labels[2], hjust = "left", size = 5)
  } else{
    p_null_neg_avg <- p_null_neg_avg +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_neg)),
                          x = uval_labels[1], y = uval_labels[2], hjust = "left", size = 5)
  }

  p_null_pos_avg <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_avg_pos_new"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Obs. - Null", y = "density", title = "Average (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_avgpos <- ggplot2::ggplot_build(p_null_pos_avg)$data[[1]]
  d_avgpos_sub <- subset(d_avgpos, .data$x > 0.1)
  if(nrow(d_avgpos_sub) > 0) {
    p_null_pos_avg <- p_null_pos_avg +
      ggplot2::geom_area(data = d_avgpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_pos_new)),
                          x = uval_labels[3], y = uval_labels[4], hjust = "left", size = 5)
  } else {
    p_null_pos_avg <- p_null_pos_avg +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_avg_pos_new)),
                          x = uval_labels[3], y = uval_labels[4], hjust = "left", size = 5)
  }
  ## Maximum
  p_null_neg_max <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_max_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Obs. - Null", y = "density", title = "Maximum (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_maxneg <- ggplot2::ggplot_build(p_null_neg_max)$data[[1]]
  d_maxneg_sub <- subset(d_maxneg, .data$x > 0.1)
  if(nrow(d_maxneg_sub) > 0) {
    p_null_neg_max <- p_null_neg_max +
      ggplot2::geom_area(data = d_maxneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_neg)),
                          x = uval_labels[5], y = uval_labels[6], hjust = "left", size = 5)
  } else {
    p_null_neg_max <- p_null_neg_max +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_neg)),
                          x = uval_labels[5], y = uval_labels[6], hjust = "left", size = 5)
  }

  p_null_pos_max <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_max_pos_new"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Obs. - Null", y = "density", title = "Maximum (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_maxpos <- ggplot2::ggplot_build(p_null_pos_max)$data[[1]]
  d_maxpos_sub <- subset(d_maxpos, .data$x > 0.1)
  if(nrow(d_maxpos_sub) > 0) {
    p_null_pos_max <- p_null_pos_max +
      ggplot2::geom_area(data = d_maxpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_pos_new)),
                          x = uval_labels[7], y = uval_labels[8], hjust = "left", size = 5)
  } else {
    p_null_pos_max <- p_null_pos_max +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_max_pos_new)),
                          x = uval_labels[7], y = uval_labels[8], hjust = "left", size = 5)
  }
  ## Variational
  p_null_neg_var <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_var_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Obs. - Null", y = "density", title = "Variational (negative)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_varneg <- ggplot2::ggplot_build(p_null_neg_var)$data[[1]]
  d_varneg_sub <- subset(d_varneg, .data$x > 0.1)
  if(nrow(d_varneg_sub) > 0) {
    p_null_neg_var <- p_null_neg_var +
      ggplot2::geom_area(data = d_varneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                 mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_var_neg)),
                 x = uval_labels[9], y = uval_labels[10], hjust = "left", size = 5)
  } else {
    p_null_neg_var <- p_null_neg_var +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_var_neg)),
                          x = uval_labels[9], y = uval_labels[10], hjust = "left", size = 5)
  }

  p_null_pos_var <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "cdelta_var_pos_new"), ggplot2::aes(x = .data$obs_minus_null)) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
    ggplot2::labs(x = "Measured unfairness", y = "density", title = "Variational (positive)") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
          legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
          axis.title.x = ggplot2::element_text(size = 12))

  d_varpos <- ggplot2::ggplot_build(p_null_pos_var)$data[[1]]
  d_varpos_sub <- subset(d_varpos, .data$x > 0.1)
  if(nrow(d_varpos_sub) > 0) {
    p_null_pos_var <- p_null_pos_var +
      ggplot2::geom_area(data = d_varpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
      ggplot2::geom_label(data = table_uval,
                 mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_var_pos_new)),
                 x = uval_labels[11], y = uval_labels[12], hjust = "left", size = 5)
  } else {
    p_null_pos_var <- p_null_pos_var +
      ggplot2::geom_label(data = table_uval,
                          mapping = ggplot2::aes(label = paste0("u-value = ", .data$cdelta_var_pos_new)),
                          x = uval_labels[11], y = uval_labels[12], hjust = "left", size = 5)
  }

  ## All in single plot
  p_null_all <- egg::ggarrange(plots = list(p_null_neg_avg, p_null_pos_avg, p_null_neg_max,
                                            p_null_pos_max, p_null_neg_var, p_null_pos_var),
                               nrow = 3, draw = F)

  # Return plots
  return(list(cfpr = p_cfpr, cfnr = p_cfnr, metrics_pos = p_metrics_pos,
              metrics_neg = p_metrics_neg, null_dist = p_null_all))
}
