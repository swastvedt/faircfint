test_that("external model of correct type", {
  pa_model_ext_nnet <- get_pa_ext_small(ex_data_external, fit_method_ext = "neural_net", maxit = 100,
                                        pa_xvars_ext = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5"))

  expect_s3_class(pa_model_ext_nnet, "nnet")

  pa_model_ext_mult <- get_pa_ext_small(ex_data_external, fit_method_ext = "multinomial", pa_xvars_ext = c("X_pa.1", "X_pa.2", "X_pa.3", "X_pa.4", "X_pa.5"))

  expect_s3_class(pa_model_ext_mult, "multinom")
})
