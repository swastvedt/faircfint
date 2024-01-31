# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })

# For formal tests: generate data with artificially set cf error rates, check defs are within a certain range

pi_model_ex <- glm(as.factor(as.character(D)) ~ A1*A2 + A1 + A2 + S_prob + X.1 + X.2 + X.3 + X.4,
                   data = ex_data_estimation, family = "binomial")

pi_model_ex_small <- glm(as.factor(as.character(D)) ~ A1*A2 + A1 + A2 + S_prob + X.1 + X.2 + X.3 + X.4,
                   data = ex_data_small, family = "binomial")

pa_model_ext_ex <- get_pa_ext_small(ex_data_external, fit_method_ext = "multinomial", pa_xvars_ext = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5"))



results_normal <- analysis_estimation(ex_data_estimation, cutoff = 0.5, estimator_type = "normal",
                    gen_null = T, bootstrap = 'rescaled', m_factor = 0.75,
                    pi_model = pi_model_ex, pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"))

results_plots <- get_plots(results_normal, sampsize = 5000, alpha = 0.05, m_factor = 0.75,
                           plot_labels = c("Group 1", "Group 2", "Group 3", "Group 4"),
                           plot_values = c(15,16,17,18), plot_colors = c("black", "gray", "red", "dodgerblue"),
                           delta_uval = 0.1)

# FIX change R and B back to 1000 when done testing
results_small <- analysis_estimation(ex_data_small, cutoff = 0.5, estimator_type = "small_internal",
                                     gen_null = T, bootstrap = 'rescaled',
                                     pi_model = pi_model_ex_small, pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"),
                                     outcome_model_type = "glm", outcome_xvars = c("X_outcome.1", "X_outcome.2", "X_outcome.3", "X_outcome.4", "X_outcome.5", "X_outcome.6", "X_outcome.7", "X_outcome.8"),
                                     fit_method_int = "multinomial", nfolds = 5, pa_xvars_int = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5", "X_pa.6", "X_pa.7"))
save(results_small, file = "C:/Users/solve/Documents/Dissertation Research/Small subgroups/rda/R package testing/results_small_all.rda")

results_small_ext <- analysis_estimation(ex_data_small, cutoff = 0.5, estimator_type = "small_borrow",
                                     gen_null = F, bootstrap = 'none',
                                     pi_model = pi_model_ex_small, pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"),
                                     outcome_model_type = "glm", outcome_xvars = c("X_outcome.1", "X_outcome.2", "X_outcome.3", "X_outcome.4", "X_outcome.5", "X_outcome.6", "X_outcome.7", "X_outcome.8"),
                                     fit_method_int = "multinomial", nfolds = 5, pa_xvars_int = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5", "X_pa.6", "X_pa.7"),
                                     data_external = ex_data_external, fit_method_ext = "multinomial", pa_xvars_ext = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5"), borrow_metric = "brier",
                                     pa_model_ext = pa_model_ext_ex)
save(results_small_ext, file = "C:/Users/solve/Documents/Dissertation Research/Small subgroups/rda/R package testing/results_small_ext_all.rda")


