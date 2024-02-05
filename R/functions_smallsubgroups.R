# Optimize multi-class AUC or Brier score for data borrowing parameter.
#
# @param alpha Borrowing parameter.
# @param preds_ext Predicted P(A=a) probabilities on internal data from external model.
# @param preds_int Predicted P(A=a) probabilities on internal data from internal model.
# @param data Internal data set.
# @param borrow_metric Which performance metric to use: Multi-class AUC ("auc") or Brier score ("brier")

# @returns Multi-class AUC or Brier score for a given value of alpha.
#
borrow_alpha <- function(alpha, preds_ext, preds_int, data, borrow_metric) {
  # Weight predicted probabilities
  preds_ext_alpha <- alpha*preds_ext
  preds_int_alpha <- (1-alpha)*preds_int
  preds_borrow <- preds_ext_alpha + preds_int_alpha

  if(borrow_metric=="auc") {
    # Multi-class AUC
    auc_val <- pROC::multiclass.roc(data$A1A2, preds_borrow)
    auc_return <- -as.numeric(auc_val$auc)
    return(auc_return)
  } else if(borrow_metric=="brier") {
    # Brier score
    brier_return <- mlr3measures::mbrier(data$A1A2, as.matrix(preds_borrow))
    return(brier_return)
  }
}

# Fit outcome and P(A=a) models, with or without data borrowing
#
# @inheritParams analysis_estimation
#
# @returns List with the following components:
#  * preds_out: List of outcome model predictions.
#  * preds_pa: Vector of P(A=a) predictions (internal-only or borrowed depending on estimator_type).
#  * alpha: Data borrowing final alpha value (returned only if estimator_type is small_borrow).
#
get_models_small <- function(data, cutoff, outcome_xvars, pa_xvars_int,
                             estimator_type, outcome_model_type, fit_method_int, nfolds,
                             data_external, fit_method_ext, borrow_metric, pa_xvars_ext,
                             pa_model_ext) {

  # Set S using cutoff
  data$S <- dplyr::if_else(data$S_prob >= cutoff, 1, 0)

  # Formulas for outcome models
  f_mu0_model <- stats::as.formula(paste0("as.factor(as.character(Y)) ~ D +", paste(outcome_xvars, collapse = " + ")))
  if(nlevels(data$A2) > 1) {
    f_mu0_A_model <- stats::as.formula(paste0("as.factor(as.character(Y)) ~ D + A1+A2+A1*A2+", paste(outcome_xvars, collapse = " + ")))
  } else {
    f_mu0_A_model <- stats::as.formula(paste0("Y ~ D + A1 +", paste(outcome_xvars, collapse = " + ")))
  }
  f_mu0_S_model <- stats::as.formula(paste0("as.factor(as.character(Y)) ~ S + D + ", paste(outcome_xvars, collapse = " + ")))
  if(nlevels(data$A2)>1) {
    f_mu0_S_A_model <- stats::as.formula(paste0("as.factor(as.character(Y)) ~ S + D + A1+A2+A1*A2 +", paste(outcome_xvars, collapse = " + ")))
  } else {
    f_mu0_S_A_model <- stats::as.formula(paste0("Y ~ S + D + A1 +", paste(outcome_xvars, collapse = " + ")))
  }
  outcome_xvars_numeric <- names(dplyr::select(dplyr::select(data, tidyr::all_of(outcome_xvars)), dplyr::where(is.numeric)))
  ## Folds for CV
  data_cv <- data
  data_cv <- dplyr::mutate(data_cv,
                    dplyr::across(tidyr::all_of(outcome_xvars_numeric), ~c(scale(.))),
                    row_id = seq(1, nrow(data_cv)),
                    A1A2 = as.factor(paste0(data$A1, data$A2)))
  data_cv <- dplyr::group_by(data_cv, .data$A1A2)
  data_cv <- dplyr::slice_sample(data_cv, prop=1)
  data_cv <- dplyr::mutate(data_cv, fold=rep(1:nfolds, length.out=dplyr::n()))
  data_cv <- dplyr::ungroup(data_cv)
  data_cv <- dplyr::arrange(data_cv, .data$row_id)
  data_cv <- dplyr::mutate(data_cv, -.data$row_id)

  # Build outcome models according to specified type
  if(outcome_model_type == "glm") {
    ### Outcome model, no S
    mu0_model <- stats::glm(f_mu0_model, data = data, family = "binomial")
    mu0hat <- prob_trunc(stats::predict(mu0_model, newdata = dplyr::mutate(data, D = 0), type = "response"))
    ### Outcome model, including S
    mu0_S_model <- stats::glm(f_mu0_S_model, data = data, family = "binomial")
    mu0S1hat <- prob_trunc(stats::predict(mu0_S_model, newdata = dplyr::mutate(data, D = 0, S = 1), type = "response"))
    mu0S0hat <- prob_trunc(stats::predict(mu0_S_model, newdata = dplyr::mutate(data, D = 0, S = 0), type = "response"))

  } else if(outcome_model_type == "neural_net") {
    mu0hat <- mu0S1hat <- mu0S0hat <- vector(mode = "numeric", length = nrow(data_cv))
    ## Fit and predict using CV
    for(i in 1:nfolds) {
      testIndexes <- which(data_cv$fold==i,arr.ind=TRUE)
      testData <- data_cv[testIndexes, ]
      trainData <- data_cv[-testIndexes, ]
      ### Outcome model, no S
      mu0_model <- nnet::nnet(f_mu0_model, data = trainData,
                              size = 100,
                              decay = 1,
                              maxit = 3000,
                              MaxNWts = 6000, trace = F)
      if(mu0_model$convergence==0) {message("Outcome model, mu0_star: converged")}else{message("Outcome model, mu0_star: did not converge")}

      mu0hat[testIndexes] <- prob_trunc(stats::predict(mu0_model, newdata = dplyr::mutate(testData, D = 0),
                                                type = "raw")[,1])
      ### Outcome model, including S
      mu0_S_model <- nnet::nnet(f_mu0_S_model, data = trainData,
                                size = 100,
                                decay = 1,
                                maxit = 3000,
                                MaxNWts = 6000, trace = F)
      if(mu0_S_model$convergence==0) {message("Outcome model, mu0: converged")}else{message("Outcome model, mu0: did not converge")}

      mu0S1hat[testIndexes] <- prob_trunc(stats::predict(mu0_S_model, newdata = dplyr::mutate(testData, D = 0, S=1),
                                                  type = "raw")[,1])
      mu0S0hat[testIndexes] <- prob_trunc(stats::predict(mu0_S_model, newdata = dplyr::mutate(testData, D = 0, S=0),
                                                  type = "raw")[,1])
    }
  }
  pred_out_list <- list("mu0hat" = mu0hat, "mu0S1hat" = mu0S1hat, "mu0S0hat" = mu0S0hat)

  # Get internal P(A=a) probabilities
  pahat_mat_int <- get_pa_int_small(data = data, pa_xvars_int = pa_xvars_int, fit_method_int = fit_method_int, nfolds = nfolds)
  ## Truncate predicted probabilities
  pahat_mat_int[] <- prob_trunc(pahat_mat_int)

  if(estimator_type == "small_borrow") {
    # Get external P(A=a) probabilities
    if(is.null(pa_model_ext)) {
      pa_model_ext <- get_pa_ext_small(data_external = data_external, pa_xvars_ext = pa_xvars_ext, fit_method_ext = fit_method_ext)
    }
    if (fit_method_ext == "multinomial") {
      pahat_mat_ext <- stats::predict(pa_model_ext, newdata = data, type = "prob")

    } else if (fit_method_ext == "neural_net") {
      pahat_mat_ext <- stats::predict(pa_model_ext, newdata = data, type = "raw")
    }
    ## Truncate
    pahat_mat_ext[] <- prob_trunc(pahat_mat_ext)

    ## Data borrowing
    data$A1A2 <- as.factor(paste0(data$A1, data$A2))
    result_temp <- stats::optim(par = 0.5, fn = borrow_alpha, method = "Brent",
                         lower = 0, upper = 1,
                         preds_ext = pahat_mat_ext, preds_int = pahat_mat_int,
                         data = data, borrow_metric = borrow_metric)
    alpha_temp <- result_temp$par

    pahat_mat <- as.matrix(as.data.frame(alpha_temp*pahat_mat_ext + (1-alpha_temp)*pahat_mat_int))

    return_list <- list(preds_out = pred_out_list, preds_pa = pahat_mat, alpha = alpha_temp)

  } else if(estimator_type == "small_internal") {
    # If no borrowing, use P(A=a) from internal data model
    pahat_mat <- pahat_mat_int
    return_list <- list(preds_out = pred_out_list, preds_pa = pahat_mat)
  }
  return(return_list)
}

# Fit and predict on internal data P(A=a) model.
#
# @inheritParams analysis_estimation
#
# @returns Predicted P(A=a|X) matrix.
#
get_pa_int_small <- function(data, pa_xvars_int, fit_method_int, nfolds) {

  pa_xvars_int_numeric <- names(dplyr::select(dplyr::select(data, tidyr::all_of(pa_xvars_int)), dplyr::where(is.numeric)))

  ## Folds for CV
  data_cv <- data
  data_cv <- dplyr::mutate(data_cv, dplyr::across(tidyr::all_of(pa_xvars_int_numeric), ~c(scale(.))))
  data_cv <- dplyr::mutate(data_cv, row_id = seq(1, nrow(data_cv)),
                           A1A2 = as.factor(paste0(data$A1, data$A2)))

  data_cv <- dplyr::group_by(data_cv, .data$A1A2)
  data_cv <- dplyr::slice_sample(data_cv, prop=1)
  data_cv <- dplyr::mutate(data_cv, fold=rep(1:nfolds, length.out=dplyr::n()))
  data_cv <- dplyr::ungroup(data_cv)
  data_cv <- dplyr::arrange(data_cv, .data$row_id)
  data_cv <- dplyr::mutate(data_cv, -.data$row_id)

  ## Fit and predict using CV
  pa_pred_int <- matrix(nrow = nrow(data_cv), ncol = nlevels(data_cv$A1A2))
  f_pa_model_int <- stats::as.formula(paste0("as.factor(as.character(A1A2)) ~ ", paste(pa_xvars_int, collapse = " + ")))

  if (fit_method_int == "multinomial") {
    for(i in 1:nfolds) {
      testIndexes <- which(data_cv$fold==i,arr.ind=TRUE)
      testData <- data_cv[testIndexes, ]
      trainData <- data_cv[-testIndexes, ]
      # (First level of A1 is baseline)
      if(nlevels(data$A1)==2 & nlevels(data$A2)==1) {
        pa_model_int <- stats::glm(f_pa_model_int, data = trainData, family = "binomial")
        pa_pred_int[testIndexes,2] <- stats::predict(pa_model_int, newdata = testData, type = "response")
        pa_pred_int[testIndexes,1] <- 1-pa_pred_int[testIndexes,2]
        colnames(pa_pred_int) <- levels(data$A1A2)
      } else {
        pa_model_int <- nnet::multinom(f_pa_model_int, data = trainData, trace = F)
        pa_pred_int[testIndexes,] <- stats::predict(pa_model_int, newdata = testData, "probs")
        colnames(pa_pred_int) <- levels(data$A1A2)
      }
    }
  } else if (fit_method_int == "neural_net") {
    for(i in 1:nfolds) {
      testIndexes <- which(data_cv$fold==i,arr.ind=TRUE)
      testData <- data_cv[testIndexes, ]
      trainData <- data_cv[-testIndexes, ]
      pa_model_int <- nnet::nnet(f_pa_model_int, data = trainData,
                                 size = 50,
                                 decay = 1,
                                 maxit = 3000,
                                 MaxNWts = 6000, trace = F)

      if(nlevels(data$A1)==2 & nlevels(data$A2)==1) {
        pa_pred_int[testIndexes,2] <- stats::predict(pa_model_int, newdata = testData,
                                              type = "raw")
        pa_pred_int[testIndexes,1] <- 1-pa_pred_int[testIndexes,2]
      } else {
        pa_pred_int[testIndexes,] <- stats::predict(pa_model_int, newdata = testData,
                                             type = "raw")
      }
      if(pa_model_int$convergence==0) {message("Internal model: converged")}else{message("Internal model: did not converge")}
    }
    colnames(pa_pred_int) <- levels(data$A1A2)
  }
  return(pa_pred_int)
}

#' Fit external data P(A=a) model.
#'
#' @inheritParams analysis_estimation
#' @param maxit Maximum number of iterations for neural_net model. Default 1000.
#'
#' @returns Model for P(A=a|X) fit on external data.
#'
#' @export

get_pa_ext_small <- function(data_external, pa_xvars_ext, fit_method_ext, maxit = 1000) {

  f_pa_model_ext <- stats::as.formula(paste0("as.factor(as.character(A1A2)) ~ ", paste(pa_xvars_ext, collapse = " + ")))
  data_external <- dplyr::mutate(data_external,
                     A1A2 = as.factor(paste0(data_external$A1, data_external$A2)))
  pa_xvars_ext_numeric <- names(dplyr::select(dplyr::select(data_external, tidyr::all_of(pa_xvars_ext)), dplyr::where(is.numeric)))

  # Fit model (no CV needed because we will predict on internal data)
  if (fit_method_ext == "multinomial") {
    pa_model_ext <- nnet::multinom(f_pa_model_ext, data = dplyr::mutate(data_external, dplyr::across(tidyr::all_of(pa_xvars_ext_numeric), ~c(scale(.)))),
                                   family = 'multinomial')
  } else if (fit_method_ext == "neural_net") {
    pa_model_ext <- nnet::nnet(f_pa_model_ext, data = dplyr::mutate(data_external, dplyr::across(tidyr::all_of(pa_xvars_ext_numeric), ~c(scale(.)))),
                               size = 100,
                               decay = 1,
                               maxit = maxit,
                               MaxNWts = 6000, trace=F)
    if(pa_model_ext$convergence==0) {message("External model: converged")}else{message("External model: did not converge")}
  }
  return(pa_model_ext)
}

