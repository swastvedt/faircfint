#' Plot results of 'analysis_estimation'.
#'
#' @param results Results of 'analysis_estimation', including null distribution and bootstrap.
#'  Must be a list with at least components 'defs' and 'boot_out'. Optional component: 'table_null'.
#' @param sampsize Sample size of estimation data set.
#' @param alpha Size of confidence interval.
#' @param m_factor Fractional power for calculating resample size.
#' @param plot_labels Character vector of protected group names for labeling counterfactual error
#'  rate plots. Must have length(A1)*length(A2) components. Default is "Group X" with groups in factor level order, varying first A2 then A1.
#' @param plot_values Numeric vector of ggplot2 shape scale values
#'  for counterfactual error rate plots. Must have length(A1)*length(A2) components. Default is all circles (shape 16).
#' @param plot_colors Character vector of colors (string name or hex) for counterfactual error
#'  rate plots. Must have length(A1)*length(A2) components. Default is viridis discrete scale.
#' @param delta_uval Threshold for u-value test. Must be between 0 and 1, default 0.1.
#'
#' @returns List of plots with the following components:
#'  * cfpr: Counterfactual false positive rates by group.
#'  * cfpr_unlabeled: Counterfactual false positive rate plot without group labels.
#'  * cfnr: Counterfactual false negative rates by group.
#'  * cfnr_unlabeled: Counterfactual false negative rate plot without group labels.
#'  * metrics_pos: Positive unfairness metrics.
#'  * metrics_neg: Negative unfairness metrics.
#'  * null_dist: Null distribution plots with estimated metrics and confidence intervals. Returned only if 'table_null' is in 'results'.
#'
#' @export

get_plots <- function(results, sampsize, alpha, m_factor,
                      plot_labels=NULL, plot_values=NULL, plot_colors=NULL, delta_uval=0.1) {
  # Prepare results
  ## Check inputs
  if('table_null' %in% names(results) & is.null(delta_uval)) stop("Must specify 'delta_uval' if null distribution results are provided.")
  est_named <- results$defs

  ## Rescaled bootstrap results
  bs_rescaled <- get_bs_rescaled(bs_table = results$boot_out, est_vals = results$defs,
                                 sampsize = sampsize, m_factor = m_factor)
  ## Confidence intervals
  ### For summary metrics
  #### t-interval
  ci_t_avg_neg <- ci_tint(bs_rescaled, est_named, parameter = "avg_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_avg_pos <- ci_tint(bs_rescaled, est_named, parameter = "avg_pos", sampsize = sampsize, alpha = alpha, m_factor = m_factor)

  ci_t_max_neg <- ci_tint(bs_rescaled, est_named, parameter = "max_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_max_pos <- ci_tint(bs_rescaled, est_named, parameter = "max_pos", sampsize = sampsize, alpha = alpha, m_factor = m_factor)

  ci_t_var_neg <- ci_tint(bs_rescaled, est_named, parameter = "var_neg", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  ci_t_var_pos <- ci_tint(bs_rescaled, est_named, parameter = "var_pos", sampsize = sampsize, alpha = alpha, m_factor = m_factor)
  #### t-interval, truncated
  ci_trunc_avg_neg <- ci_trunc(ci_t_avg_neg, type = "tint")
  ci_trunc_avg_pos <- ci_trunc(ci_t_avg_pos, type = "tint")

  ci_trunc_max_neg <- ci_trunc(ci_t_max_neg, type = "tint")
  ci_trunc_max_pos <- ci_trunc(ci_t_max_pos, type = "tint")

  ci_trunc_var_neg <- ci_trunc(ci_t_var_neg, type = "tint")
  ci_trunc_var_pos <- ci_trunc(ci_t_var_pos, type = "tint")

  ci_trunc_df <- list("avg_neg" = ci_trunc_avg_neg, "avg_pos" = ci_trunc_avg_pos,
                      "max_neg" = ci_trunc_max_neg, "max_pos" = ci_trunc_max_pos,
                      "var_neg" = ci_trunc_var_neg, "var_pos" = ci_trunc_var_pos) %>% dplyr::bind_rows(.id = "stat")
  rownames(ci_trunc_df) <- NULL
  ### For counterfactual error rates
  results_cfpr_temp <- est_named[grep("cfpr", names(est_named))]
  results_cfpr <- results_cfpr_temp[-grep("marg", names(results_cfpr_temp))]

  results_cfnr_temp <- est_named[grep("cfnr", names(est_named))]
  results_cfnr <- results_cfnr_temp[-grep("marg", names(results_cfnr_temp))]

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

  ## Named null distribution results
  if('table_null' %in% names(results)) {
    table_null <- results$table_null

    est_named_pivot <- data.frame(t(est_named)) %>% dplyr::select(.data$avg_neg, .data$avg_pos,
                                                   .data$max_neg, .data$max_pos,
                                                   .data$var_neg, .data$var_pos) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "stat", values_to = "value_obs")

    table_null_delta <- as.data.frame(table_null) %>% dplyr::select(.data$avg_neg, .data$avg_pos,
                                                     .data$max_neg, .data$max_pos,
                                                     .data$var_neg, .data$var_pos) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "stat", values_to = "value_null") %>%
      dplyr::left_join(est_named_pivot, by = "stat") %>%
      dplyr::mutate(obs_minus_null = .data$value_obs - .data$value_null)

    ## U-value table
    table_uval <- data.frame(
      avg_neg = mean((est_named["avg_neg"] - table_null[,"avg_neg"]) > delta_uval, na.rm = T),
      avg_pos = mean((est_named["avg_pos"] - table_null[,"avg_pos"]) > delta_uval, na.rm = T),
      max_neg = mean((est_named["max_neg"] - table_null[,"max_neg"]) > delta_uval, na.rm = T),
      max_pos = mean((est_named["max_pos"] - table_null[,"max_pos"]) > delta_uval, na.rm = T),
      var_neg = mean((est_named["var_neg"] - table_null[,"var_neg"]) > delta_uval, na.rm = T),
      var_pos = mean((est_named["var_pos"] - table_null[,"var_pos"]) > delta_uval, na.rm = T)
    ) %>% cond_round_3()
  }

  ## Estimation table with confidence intervals
  ### Put all CI information in same table
  allcis_trunc <- dplyr::bind_rows(list(ci_trunc_df, ci_trunc_cfpr, ci_trunc_cfnr))
  rownames(allcis_trunc) <- NULL
  ### Add CIs to estimation table
  est_summaries <- data.frame(t(est_named)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "stat") %>%
    dplyr::mutate(sign = dplyr::case_when(
      stringr::str_detect(.data$stat, "cfpr_marg") ~ "cfpr_marg",
      stringr::str_detect(.data$stat, "cfnr_marg") ~ "cfnr_marg",
      stringr::str_detect(.data$stat, "cfpr") ~ "cfpr",
      stringr::str_detect(.data$stat, "cfnr") ~ "cfnr",
      stringr::str_detect(.data$stat, "fpr") ~ "fpr",
      stringr::str_detect(.data$stat, "fnr") ~ "fnr",
      .data$stat %in% c("avg_neg", "max_neg", "var_neg") ~ "aggregate_neg",
      .data$stat %in% c("avg_pos", "max_pos", "var_pos") ~ "aggregate_pos",
      TRUE ~ "other"
    )) %>% dplyr::left_join(allcis_trunc, by = "stat")

  # Plots
  ## Set defaults if not user-specified
  num_gps <- length(grep("cfpr", names(results$defs))) - length(grep("cfpr_marg", names(results$defs)))

  if(is.null(plot_labels)) {
    gp_names <- stringr::str_split(names(results$defs)[1:num_gps], "_", simplify = T)[,2]
    plot_labels <- paste0("Group ", gp_names)
  }
  if(is.null(plot_values)) { plot_values <- rep(16, num_gps) }

  ## Counterfactual error rate plots
  p_cfpr <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfpr"), ggplot2::aes(x = "cFPR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFPR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_shape_manual(values = plot_values, labels = plot_labels)

  p_cfpr_unlabeled <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfpr"), ggplot2::aes(x = "cFPR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFPR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_shape_manual(values = plot_values)

  p_cfnr <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfnr"), ggplot2::aes(x = "cFNR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFNR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
          legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_shape_manual(values = plot_values, labels = plot_labels)

  p_cfnr_unlabeled <- ggplot2::ggplot(dplyr::filter(est_summaries, .data$sign == "cfnr"), ggplot2::aes(x = "cFNR", y = .data$value, color = .data$stat, shape = .data$stat)) +
    ggplot2::geom_point(size = 4, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$low_trans, ymax = .data$high_trans), width = 0.2, position=ggplot2::position_dodge(width=0.4)) +
    ggplot2::labs(x = NULL, y = "Group cFNR Estimate", title = NULL, shape = "Group", color = "Group") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 14),
                   legend.text = ggplot2::element_text(size = 14), legend.spacing.y = ggplot2::unit(0.3, 'cm')) +
    ggplot2::guides(shape = ggplot2::guide_legend(byrow = TRUE)) +
    ggplot2::scale_shape_manual(values = plot_values)

  ### Set colors if not user-specified
  if(is.null(plot_colors)) {
    p_cfpr <- p_cfpr + ggplot2::scale_color_viridis_d(labels = plot_labels)
    p_cfpr_unlabeled <- p_cfpr_unlabeled + ggplot2::scale_color_viridis_d()
    p_cfnr <- p_cfnr + ggplot2::scale_color_viridis_d(labels = plot_labels)
    p_cfnr_unlabeled <- p_cfnr_unlabeled + ggplot2::scale_color_viridis_d()
  } else {
    p_cfpr <- p_cfpr + ggplot2::scale_color_discrete(type = plot_colors, labels = plot_labels)
    p_cfpr_unlabeled <- p_cfpr_unlabeled + ggplot2::scale_color_discrete(type = plot_colors)
    p_cfnr <- p_cfnr + ggplot2::scale_color_discrete(type = plot_colors, labels = plot_labels)
    p_cfnr_unlabeled <- p_cfnr_unlabeled + ggplot2::scale_color_discrete(type = plot_colors)
  }
  ## Unfairness metric plots
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

  ## Null distribution/estimation plots
  if('table_null' %in% names(results)) {
    ### Average
    p_null_neg_avg <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "avg_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Obs. - Null", y = "density", title = "Average (negative)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_avgneg <- ggplot2::ggplot_build(p_null_neg_avg)$data[[1]]
    d_avgneg_sub <- dplyr::filter(d_avgneg, .data$x > 0.1)
    uval_avgneg_x <- max(d_avgneg$x)
    uval_avgneg_y <- stats::quantile(d_avgneg$y, 0.25)

    if(nrow(d_avgneg_sub) > 0) {
      p_null_neg_avg <- p_null_neg_avg +
        ggplot2::geom_area(data = d_avgneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$avg_neg)),
                            x = uval_avgneg_x, y = uval_avgneg_y, hjust = "right", vjust = "bottom", size = 5)
    } else{
      p_null_neg_avg <- p_null_neg_avg +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$avg_neg)),
                            x = uval_avgneg_x, y = uval_avgneg_y, hjust = "right", vjust = "bottom", size = 5)
    }

    p_null_pos_avg <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "avg_pos"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Obs. - Null", y = "density", title = "Average (positive)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_avgpos <- ggplot2::ggplot_build(p_null_pos_avg)$data[[1]]
    d_avgpos_sub <- dplyr::filter(d_avgpos, .data$x > 0.1)
    uval_avgpos_x <- max(d_avgpos$x)
    uval_avgpos_y <- stats::quantile(d_avgpos$y, 0.25)

    if(nrow(d_avgpos_sub) > 0) {
      p_null_pos_avg <- p_null_pos_avg +
        ggplot2::geom_area(data = d_avgpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$avg_pos)),
                            x = uval_avgpos_x, y = uval_avgpos_y, hjust = "right", vjust = "bottom", size = 5)
    } else {
      p_null_pos_avg <- p_null_pos_avg +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$avg_pos)),
                            x = uval_avgpos_x, y = uval_avgpos_y, hjust = "right", vjust = "bottom", size = 5)
    }
    ### Maximum
    p_null_neg_max <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "max_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Obs. - Null", y = "density", title = "Maximum (negative)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_maxneg <- ggplot2::ggplot_build(p_null_neg_max)$data[[1]]
    d_maxneg_sub <- dplyr::filter(d_maxneg, .data$x > 0.1)
    uval_maxneg_x <- max(d_maxneg$x)
    uval_maxneg_y <- stats::quantile(d_maxneg$y, 0.25)

    if(nrow(d_maxneg_sub) > 0) {
      p_null_neg_max <- p_null_neg_max +
        ggplot2::geom_area(data = d_maxneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$max_neg)),
                            x = uval_maxneg_x, y = uval_maxneg_y, hjust = "right", vjust = "bottom", size = 5)
    } else {
      p_null_neg_max <- p_null_neg_max +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$max_neg)),
                            x = uval_maxneg_x, y = uval_maxneg_y, hjust = "right", vjust = "bottom", size = 5)
    }

    p_null_pos_max <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "max_pos"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Obs. - Null", y = "density", title = "Maximum (positive)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_maxpos <- ggplot2::ggplot_build(p_null_pos_max)$data[[1]]
    d_maxpos_sub <- dplyr::filter(d_maxpos, .data$x > 0.1)
    uval_maxpos_x <- max(d_maxpos$x)
    uval_maxpos_y <- stats::quantile(d_maxpos$y, 0.25)

    if(nrow(d_maxpos_sub) > 0) {
      p_null_pos_max <- p_null_pos_max +
        ggplot2::geom_area(data = d_maxpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$max_pos)),
                            x = uval_maxpos_x, y = uval_maxpos_y, hjust = "right", vjust = "bottom", size = 5)
    } else {
      p_null_pos_max <- p_null_pos_max +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$max_pos)),
                            x = uval_maxpos_x, y = uval_maxpos_y, hjust = "right", vjust = "bottom", size = 5)
    }
    ### Variational
    p_null_neg_var <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "var_neg"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Obs. - Null", y = "density", title = "Variational (negative)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_varneg <- ggplot2::ggplot_build(p_null_neg_var)$data[[1]]
    d_varneg_sub <- dplyr::filter(d_varneg, .data$x > 0.1)
    uval_varneg_x <- max(d_varneg$x)
    uval_varneg_y <- stats::quantile(d_varneg$y, 0.25)

    if(nrow(d_varneg_sub) > 0) {
      p_null_neg_var <- p_null_neg_var +
        ggplot2::geom_area(data = d_varneg_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$var_neg)),
                            x = uval_varneg_x, y = uval_varneg_y, hjust = "right", vjust = "bottom", size = 5)
    } else {
      p_null_neg_var <- p_null_neg_var +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$var_neg)),
                            x = uval_varneg_x, y = uval_varneg_y, hjust = "right", vjust = "bottom", size = 5)
    }

    p_null_pos_var <- ggplot2::ggplot(dplyr::filter(table_null_delta, .data$stat == "var_pos"), ggplot2::aes(x = .data$obs_minus_null)) +
      ggplot2::geom_density() +
      ggplot2::geom_vline(xintercept = delta_uval, linetype = 2) +
      ggplot2::labs(x = "Measured unfairness", y = "density", title = "Variational (positive)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 11), axis.text.x = ggplot2::element_text(size = 9),
                     axis.title.x = ggplot2::element_text(size = 12))

    d_varpos <- ggplot2::ggplot_build(p_null_pos_var)$data[[1]]
    d_varpos_sub <- dplyr::filter(d_varpos, .data$x > 0.1)
    uval_varpos_x <- max(d_varpos$x)
    uval_varpos_y <- stats::quantile(d_varpos$y, 0.25)

    if(nrow(d_varpos_sub) > 0) {
      p_null_pos_var <- p_null_pos_var +
        ggplot2::geom_area(data = d_varpos_sub, ggplot2::aes(x=.data$x, y=.data$y), fill = "#d3d3d3") +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$var_pos)),
                            x = uval_varpos_x, y = uval_varpos_y, hjust = "right", vjust = "bottom", size = 5)
    } else {
      p_null_pos_var <- p_null_pos_var +
        ggplot2::geom_label(data = table_uval,
                            mapping = ggplot2::aes(label = paste0("u-value = ", .data$var_pos)),
                            x = uval_varpos_x, y = uval_varpos_y, hjust = "right", vjust = "bottom", size = 5)
    }
    ### All in single plot
    p_null_all <- egg::ggarrange(plots = list(p_null_neg_avg, p_null_pos_avg, p_null_neg_max,
                                              p_null_pos_max, p_null_neg_var, p_null_pos_var),
                                 nrow = 3, draw = F)
  }
  # Return plots
  return_list <- list(cfpr = p_cfpr, cfnr = p_cfnr, cfpr_unlabeled = p_cfpr_unlabeled, cfnr_unlabeled = p_cfnr_unlabeled,
                      metrics_pos = p_metrics_pos, metrics_neg = p_metrics_neg)
  if('table_null' %in% names(results)) {
    return_list[['null_dist']] <- p_null_all
  }
  return(return_list)
}
