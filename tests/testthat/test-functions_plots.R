test_that("plots returned", {
  pi_model_ex <- glm(as.factor(as.character(D)) ~ A1*A2 + A1 + A2 + S_prob + X.1 + X.2 + X.3 + X.4,
                     data = ex_data_estimation, family = "binomial")
  results_standard <- analysis_estimation(ex_data_estimation, cutoff = 0.5, estimator_type = "standard",
                                        gen_null = T, bootstrap = 'rescaled', m_factor = 0.75,
                                        pi_model = pi_model_ex, pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"))

  results_plots <- get_plots(results_standard, sampsize = 5000, alpha = 0.05, m_factor = 0.75,
                             plot_labels = c("Group 1", "Group 2", "Group 3", "Group 4"),
                             plot_values = c(15,16,17,18), plot_colors = c("black", "gray", "red", "dodgerblue"),
                             delta_uval = 0.1)

  expect_true(ggplot2::is.ggplot(results_plots$cfpr))
  expect_true(ggplot2::is.ggplot(results_plots$cfnr))
  expect_true(ggplot2::is.ggplot(results_plots$cfpr_unlabeled))
  expect_true(ggplot2::is.ggplot(results_plots$cfnr_unlabeled))
  expect_true(ggplot2::is.ggplot(results_plots$metrics_pos))
  expect_true(ggplot2::is.ggplot(results_plots$metrics_neg))
})
