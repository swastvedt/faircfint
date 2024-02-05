test_that("cf defs return correct values", {
  data_fixed <- as.data.frame(data_gen_fixedcf(50000))
  defs_fixed <- get_defs_analysis(dplyr::mutate(data_fixed, pi = 0),
                                  cutoff = 0.5, estimator_type = "standard")
  cfpr_vals <- defs_fixed$defs[9:12]
  cfnr_vals <- defs_fixed$defs[13:16]

  expect_lt(cfpr_vals[1], 0.14)
  expect_lt(cfpr_vals[2], 0.24)
  expect_lt(cfpr_vals[3], 0.19)
  expect_lt(cfpr_vals[4], 0.29)

  expect_gt(cfpr_vals[1], 0.06)
  expect_gt(cfpr_vals[2], 0.16)
  expect_gt(cfpr_vals[3], 0.11)
  expect_gt(cfpr_vals[4], 0.21)

  expect_lt(cfnr_vals[1], 0.24)
  expect_lt(cfnr_vals[2], 0.34)
  expect_lt(cfnr_vals[3], 0.29)
  expect_lt(cfnr_vals[4], 0.39)

  expect_gt(cfnr_vals[1], 0.16)
  expect_gt(cfnr_vals[2], 0.26)
  expect_gt(cfnr_vals[3], 0.21)
  expect_gt(cfnr_vals[4], 0.31)
})

test_that("error rate 0/1/NA warnings", {
  data_toosmall <- data.frame(
    A1 = as.factor(c(rep(1,10), rep(0,10))),
    A2 = as.factor(c(rep(0,5), rep(1,5), rep(0,5), rep(1,5))),
    Y = rep(c(0,1), 10),
    Y0 = rep(c(0,1), 10),
    D = rep(0,20),
    pi = rep(0,20),
    S_prob = c(rep(0,10), rep(1,10))
  )
  expect_warning(
    expect_warning(
      expect_warning(
        expect_warning(get_defs_analysis(data_toosmall, cutoff = 0.5, estimator_type = "standard"),
                       "cFPR of 0, 1, or NULL for at least one group."),
        "cFNR of 0, 1, or NULL for at least one group."
      ),
      "Observational FPR of 0, 1, or NULL for at least one group."
    ),
    "Observational FNR of 0, 1, or NULL for at least one group."
  )
})

test_that("Y and D factor level conversion", {
  data_YD <- as.data.frame(data_gen_fixedcf(1000))
  # Recode Y and D levels
  data_YD_recode <- data_YD
  data_YD_recode$Y <- dplyr::case_when(
    data_YD$Y==0 ~ 1,
    data_YD$Y==1 ~ 2,
  )
  data_YD_recode$D <- dplyr::case_when(
    data_YD$D==0 ~ 1,
    data_YD$D==1 ~ 2,
  )

  # Check function resets to 0/1
  result_1 <- analysis_estimation(data_YD, cutoff = 0.5, estimator_type = "standard",
                                  gen_null = F, bootstrap = 'none',
                                  pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"))
  result_2 <- analysis_estimation(data_YD_recode, cutoff = 0.5, estimator_type = "standard",
                                  gen_null = F, bootstrap = 'none',
                                  pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"))
  expect_equal(result_1$defs, result_2$defs)
})

test_that("missing data warning", {
  data_miss <- as.data.frame(data_gen_fixedcf(1000))

  miss_ind <- sample(1:1000, 10, replace = F)
  data_miss[miss_ind,"X.1"] <- NA_real_
  miss_ind_ext <- sample(1:nrow(ex_data_external), 5, replace = F)
  data_miss_ext <- ex_data_external
  data_miss_ext[miss_ind_ext,"A1"] <- NA_real_

  expect_warning(analysis_estimation(data_miss, cutoff = 0.5, estimator_type = "standard",
                                  gen_null = F, bootstrap = 'none',
                                  pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4")),
                 "'data' contains missing values")
  expect_warning(analysis_estimation(ex_data_small, cutoff = 0.5, estimator_type = "small_borrow",
                                     gen_null = F, bootstrap = 'none',
                                     pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"),
                                     outcome_model_type = "glm", outcome_xvars = c("X_outcome.1", "X_outcome.2", "X_outcome.3", "X_outcome.4", "X_outcome.5", "X_outcome.6", "X_outcome.7", "X_outcome.8"),
                                     fit_method_int = "multinomial", nfolds = 5, pa_xvars_int = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5", "X_pa.6", "X_pa.7"),
                                     data_external = data_miss_ext, fit_method_ext = "multinomial", pa_xvars_ext = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5"), borrow_metric = "brier"),
                 "External data contains missing values")

})

