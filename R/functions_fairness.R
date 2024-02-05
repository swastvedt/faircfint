prob_trunc <- function(p) pmax(pmin(p, 0.995),0.005)

# Get estimated unfairness metrics.
#
# 'get_defs_analysis' returns a list of unfairness definitions.
#
# Note: Indices and/or order of definitions in output of this function are
#  used in formatting functions.
#
# @param data Dataframe or tibble. Must include at least A1, A2, Y (binary 0/1),
#  D (binary 0/1), S_prob (probability), Y0 (estimated), pi (estimated).
# @param cutoff Classification cutoff. Must be a single numeric value.
# @param estimator_type Type of estimation: standard ('standard'), small subgroup internal-only ('small_internal'),
#  or small subgroup with data borrowing ('small_borrow').
# @param preds_out Outcome model predictions (required if estimator_type is 'small_internal' or 'small_borrow').
# @param preds_pa P(A=a) model predictions (required if estimator_type is 'small_internal' or 'small_borrow').
#
# @returns
# A list with the following components:
#
# * defs: vector of unfairness definitions with 12+(4*nlevels(A1)nlevels(A2))+2*nlevels(A1)+2*nlevels(A2) items

get_defs_analysis <- function(data, cutoff, estimator_type,
                              preds_out=NULL, preds_pa=NULL) {
  # Grid of values for protected characteristic combinations
  # SMALL SUBGROUPS UPDATE: Reordering grid to match ordering of the factor A1A2 (pasted A1, A2). Previously A1 varied first, then A2.
  if(nlevels(data$A2)>1) {
    A1_levels <- levels(data$A1)
    A2_levels <- levels(data$A2)
    A1A2_grid <- expand.grid(A2_levels, A1_levels)
  } else {
    A1A2_grid <- data.frame(Var2 = levels(data$A1), Var1 = 1)
  }
  # Set up A1A2 factor
  for(i in 1:nrow(A1A2_grid)) {
    varname <- paste0("A1A2_", i)
    data[varname] <- dplyr::if_else(data$A1 == A1A2_grid[i,"Var2"] & data$A2 == A1A2_grid[i,"Var1"],
                             1, 0)
  }
  data$A1A2 <- as.factor(paste0(as.character(data$A1), as.character(data$A2)))

  # Other variables for fairness definitions
  Y0 <- data$Y0
  S <- dplyr::if_else(data$S_prob >= cutoff, 1, 0)
  Y <- as.numeric(data$Y)
  D <- as.numeric(data$D)
  pi <- data$pi
  pahat_mat <- preds_pa

  ### Intersectional FPRs and FNRs (observational)
  fpr_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*S*(1-Y))/mean(data[,varname]*(1-Y))
  })

  fnr_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*(1-S)*Y)/mean(data[,varname]*Y)
  })

  # Intersectional cFPRs and cFNRs
  cfpr_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*S*(1-Y)*(1-D)/(1-pi))/mean(data[,varname]*(1-Y)*(1-D)/(1-pi))
  })

  cfnr_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*(1-S)*Y0)/mean(data[,varname]*Y0)
  })

  ### Remove groups with 0 or 1 and give a warning
  cfpr_vec_dropped <- cfpr_vec[cfpr_vec != 0 & cfpr_vec != 1 & !is.null(cfpr_vec)]
  cfnr_vec_dropped <- cfnr_vec[cfnr_vec != 0 & cfnr_vec != 1 & !is.null(cfnr_vec)]
  fpr_vec_dropped <- fpr_vec[fpr_vec != 0 & fpr_vec != 1 & !is.null(fpr_vec)]
  fnr_vec_dropped <- fnr_vec[fnr_vec != 0 & fnr_vec != 1 & !is.null(fnr_vec)]

  if(length(cfpr_vec_dropped)!=length(cfpr_vec)) {
    # print(paste0("cFPR: ", cfpr_vec))
    # print(table(A1=data$A1, A2=data$A2, Y=data$Y, S=data$S, D=data$D))
    warning("cFPR of 0, 1, or NULL for at least one group.", call. = F)
  }
  if(length(cfnr_vec_dropped)!=length(cfnr_vec)) {
    # print(paste0("cFNR: ", cfnr_vec))
    # print(table(A1=data$A1, A2=data$A2, Y=data$Y, S=data$S, D=data$D))
    warning("cFNR of 0, 1, or NULL for at least one group.", call. = F)
  }
  if(length(fpr_vec_dropped)!=length(fpr_vec)) {
    warning("Observational FPR of 0, 1, or NULL for at least one group.", call. = F)
  }
  if(length(fnr_vec_dropped)!=length(fnr_vec)) {
    warning("Observational FNR of 0, 1, or NULL for at least one group.", call. = F)
  }

  # Marginal cFPRs and cFNRs
  cfpr_marg_vec_A1 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var2))), function(i) {
    varvec <- dplyr::if_else(data$A1 == levels(data$A1)[i], 1, 0)
    mean(varvec*S*(1-Y)*(1-D)/(1-pi))/mean(varvec*(1-Y)*(1-D)/(1-pi))
  })

  cfpr_marg_vec_A2 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var1))), function(i) {
    varvec <- dplyr::if_else(data$A2 == levels(data$A2)[i], 1, 0)
    mean(varvec*S*(1-Y)*(1-D)/(1-pi))/mean(varvec*(1-Y)*(1-D)/(1-pi))
  })

  cfnr_marg_vec_A1 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var2))), function(i) {
    varvec <- dplyr::if_else(data$A1 == levels(data$A1)[i], 1, 0)
    mean(varvec*(1-S)*Y0)/mean(varvec*Y0)
  })
  cfnr_marg_vec_A2 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var1))), function(i) {
    varvec <- dplyr::if_else(data$A2 == levels(data$A2)[i], 1, 0)
    mean(varvec*(1-S)*Y0)/mean(varvec*Y0)
  })

  # Intersectional deltas
  cdeltaps <- tryCatch({
    c(stats::dist(cfpr_vec_dropped, method = "manhattan"))
  }, error=function(cond) {NULL})
  cdeltans <- tryCatch({
    c(stats::dist(cfnr_vec_dropped, method = "manhattan"))
  }, error=function(cond) {NULL})

  # Average
  if(length(cdeltans)>0) {
    cdelta_avg_neg <- mean(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_avg_neg <- NA_real_ }

  if(length(cdeltaps)>0) {
    cdelta_avg_pos <- mean(sapply(cdeltaps, abs), na.rm = T)
  } else{ cdelta_avg_pos <- NA_real_ }

  # Max
  if(length(cdeltans)>0) {
    cdelta_max_neg <- max(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_max_neg <- NA_real_ }

  if(length(cdeltaps)>0) {
    cdelta_max_pos <- max(sapply(cdeltaps, abs), na.rm = T)
  } else{ cdelta_max_pos <- NA_real_ }

  # Variational
  if(length(cdeltans)>0) {
    cdelta_var_neg <- stats::var(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_var_neg <- NA_real_ }

  if(length(cdeltaps)>0) {
    cdelta_var_pos <- stats::var(sapply(cdeltaps, abs), na.rm = T)
  } else { cdelta_var_pos <- NA_real_ }

  ############# SMALL SUBGROUP METRICS ############################
  if(estimator_type %in% c("small_internal", "small_borrow")) {
    ## Overall cFPR and cFNR
    cfpr_all <- mean(S*(1-Y)*(1-D)/(1-pi))/mean((1-Y)*(1-D)/(1-pi))
    cfnr_all <- mean(((1-S)*Y*(1-D))/(1-pi))/mean((Y*(1-D))/(1-pi))
    ## Estimate numerators
    cfprnum_vec_small <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
      varname <- paste0("A1A2_", i)
      mean((1-preds_out$mu0S1hat)*S*data[,varname])/mean((1-preds_out$mu0S1hat)*S)
    })
    cfnrnum_vec_small <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
      varname <- paste0("A1A2_", i)
      mean(preds_out$mu0S0hat*(1-S)*data[,varname])/mean(preds_out$mu0S0hat*(1-S))
    })
    ## Estimate cFPR and cFNR small subgroup versions
    # (phat_mat columns are in factor level order, so are A1A2_grid rows)
    cfpr_vec_small <- sapply(seq(1,nrow(A1A2_grid)), function(r) {
      cfpr_all*cfprnum_vec_small[r]/(mean(pahat_mat[,r]*(1-preds_out$mu0hat))/mean(1-preds_out$mu0hat))
    })
    cfnr_vec_small <- sapply(seq(1,nrow(A1A2_grid)), function(r) {
      cfnr_all*cfnrnum_vec_small[r]/(mean(pahat_mat[,r]*preds_out$mu0hat)/mean(preds_out$mu0hat))
    })
    ### Remove groups with 0 or 1 and give a warning
    cfpr_vec_small_dropped <- cfpr_vec_small[cfpr_vec_small != 0 & cfpr_vec_small != 1]
    cfnr_vec_small_dropped <- cfnr_vec_small[cfnr_vec_small != 0 & cfnr_vec_small != 1]
    if(length(cfpr_vec_small_dropped)!=length(cfpr_vec_small)) {
      warning("Small subgroup cFPR of 0 or 1 for at least one group.")
    }
    if(length(cfnr_vec_small_dropped)!=length(cfnr_vec_small)) {
      warning("Small subgroup cFNR of 0 or 1 for at least one group.")
    }

    ## Intersectional deltas
    cdeltaps_small <- tryCatch({
      c(stats::dist(cfpr_vec_small_dropped, method = "manhattan"))
    }, error=function(cond) {NULL})

    cdeltans_small <- tryCatch({
      c(stats::dist(cfnr_vec_small_dropped, method = "manhattan"))
    }, error=function(cond) {NULL})

    ## Average, Max, Variational
    if(!is.null(cdeltans_small)) {
      cdelta_avg_neg_small <- mean(sapply(cdeltans_small, abs))
      cdelta_max_neg_small <- max(sapply(cdeltans_small, abs))
      cdelta_var_neg_small <- stats::var(sapply(cdeltans_small, abs))
    } else { cdelta_avg_neg_small <- cdelta_max_neg_small <- cdelta_var_neg_small <- NA_real_ }

    if(!is.null(cdeltaps_small)) {
      cdelta_avg_pos_small <- mean(sapply(cdeltaps_small, abs))
      cdelta_max_pos_small <- max(sapply(cdeltaps_small, abs))
      cdelta_var_pos_small <- stats::var(sapply(cdeltaps_small, abs))
    } else { cdelta_avg_pos_small <- cdelta_max_pos_small <- cdelta_var_pos_small <- NA_real_ }
  }

  # Return vector of metrics
  defs <- c(fpr_vec, fnr_vec, cfpr_vec, cfnr_vec,
            cfpr_marg_vec_A1, cfpr_marg_vec_A2, cfnr_marg_vec_A1, cfnr_marg_vec_A2,
            cdelta_avg_neg, cdelta_max_neg, cdelta_var_neg,
            cdelta_avg_pos, cdelta_max_pos, cdelta_var_pos
  )
  if(estimator_type %in% c("small_internal", "small_borrow")) {
    defs <- c(defs,
              cfpr_vec_small, cfnr_vec_small,
              cdelta_avg_neg_small, cdelta_max_neg_small, cdelta_var_neg_small,
              cdelta_avg_pos_small, cdelta_max_pos_small, cdelta_var_pos_small)
  }
  return(list(defs = defs))
}

# Rescaled, stratified bootstrap for unfairness metrics.
#
# @inheritParams analysis_estimation
# @param B Number of bootstrap replications.
# @param m_factor Fractional power for calculating resample size (default 0.75).
# @param defs_names Character vector of names for get_defs_analysis output.
#
# @returns A matrix of definitions, one row for each of the B bootstrap replications.

bs_rescaled_analysis <- function(data, cutoff, B, m_factor,
                                 pi_model_type, pi_model_seed = NULL, pi_xvars,
                                 defs_names, estimator_type,
                                 outcome_model_type, outcome_xvars,
                                 fit_method_int, nfolds, pa_xvars_int, data_external, fit_method_ext, borrow_metric, pa_xvars_ext,
                                 pa_model_ext) {
  # Empty matrix for results
  boot.out <- matrix(nrow = B, ncol = length(defs_names))
  colnames(boot.out) <- defs_names
  # m size for m of n bootstrap
  n_rescaled <- floor(nrow(data)^m_factor)

  for (b in 1:B) {
    # Stratified sample conditional on proportions of (A1, A2), Y, S (binary)
    data$A1A2 <- paste0(data$A1, data$A2)
    data$S_char <- dplyr::if_else(data$S_prob >= cutoff, "1", "0")

    data_bs <- as.data.frame(splitstackshape::stratified(data, c("A1A2", "Y", "S_char"), n_rescaled/nrow(data), replace = F))
    # Get Y0 estimator
    ests <- get_est_analysis(data_bs, cutoff = cutoff, pi_model = NULL,
                             pi_model_type = pi_model_type, pi_model_seed = pi_model_seed, pi_xvars = pi_xvars)
    ## Get small subgroup unfairness metrics
    if(estimator_type %in% c("small_internal", "small_borrow")) {
      ### Get outcome and P(A=a) models
      models_temp <- get_models_small(data=data_bs, cutoff=cutoff, estimator_type=estimator_type,
                                      outcome_xvars=outcome_xvars, outcome_model_type=outcome_model_type,
                                      pa_xvars_int=pa_xvars_int, fit_method_int=fit_method_int, nfolds=nfolds,
                                      data_external=data_external, fit_method_ext=fit_method_ext, borrow_metric=borrow_metric, pa_xvars_ext=pa_xvars_ext,
                                      pa_model_ext=pa_model_ext)
      ### Get metrics
      defs_temp <- get_defs_analysis(dplyr::mutate(data_bs, Y0 = ests[['Y0est']], pi = ests[['pi']]),
                                        cutoff=cutoff, estimator_type=estimator_type,
                                        preds_out = models_temp$preds_out, preds_pa = models_temp$preds_pa)[['defs']]
    } else {
      ## Get original unfairness metrics
      defs_temp <- get_defs_analysis(dplyr::mutate(data_bs, Y0 = ests[['Y0est']], pi = ests[['pi']]),
                                     cutoff=cutoff, estimator_type=estimator_type)[['defs']]
    }

    # Add row to matrix
    boot.out[b,] <- defs_temp
  }
  return(boot.out)
}

# Get estimated propensity score and $Y^0$ weighted estimator.
#
# @inheritParams analysis_estimation
#
# @returns List with the following components:
#  * Y0est: Vector of inverse probability weighted estimates of $Y^0$.
#  * pi: Vector of estimated propensity scores.
#

get_est_analysis <- function(data, cutoff, pi_model, pi_model_type, pi_model_seed, pi_xvars) {
  # Blank list for storing results
  est_choice <- list("Y0est" = NULL, "pi" = NULL)
  # Classify predictions according to cutoff
  data$S <- dplyr::if_else(data$S_prob >= cutoff, 1, 0)

  ## Get propensity score model if pi_model not specified
  if(is.null(pi_model)) {
    if(pi_model_type == "glm") {
      pi_model_formula <- paste0("D ~ ", paste(pi_xvars[1:length(pi_xvars)-1], collapse=" + "), " + ", pi_xvars[length(pi_xvars)])
      pi_model <- stats::glm(pi_model_formula, data = data, family = "binomial")

    } else if(pi_model_type == "rf") {
      data_x <- dplyr::select(data, dplyr::all_of(pi_xvars))
      matrix_x <- stats::model.matrix(~.-1, data = data_x)

      set.seed(pi_model_seed[1])
      pi_model <- randomForest::randomForest(matrix_x, as.factor(data$D), importance = T,
                               ntree = 2000, mtry = 90, nodesize = 50)
    }
  }
  # Get propensity score predictions
  if(pi_model_type == 'glm') {
    pihat <- prob_trunc(stats::predict(pi_model, newdata = data, type = "response"))
  } else if(pi_model_type == 'rf') {
    set.seed(pi_model_seed[2])
    pihat <- prob_trunc(stats::predict(pi_model, newdata = matrix_x, type = "prob")[,2])
  } else stop("'pi_model_type' must be 'glm' or 'rf'")
  ## IPW estimate of Y0
  ipwhat <- ((1-data$D)*data$Y)/(1-pihat)
  ## Save Y0 and pi estimates
  est_choice[['Y0est']] <- ipwhat
  est_choice[['pi']] <- pihat

  return(est_choice)
}

# Simulate null distributions.
#
# @inheritParams bs_rescaled_analysis
# @param R Number of replications.
#
# @returns List with the following components:
#  * ipw: Matrix of null unfairness metrics, one row for each replication.

analysis_nulldist <- function(data, R, cutoff,
                              pi_model_type, pi_model_seed, pi_xvars,
                              defs_names, estimator_type,
                              outcome_model_type, outcome_xvars,
                              fit_method_int, nfolds, pa_xvars_int, data_external, fit_method_ext, borrow_metric, pa_xvars_ext,
                              pa_model_ext) {
  ## Matrix for storing results
  table.null <- matrix(nrow = R, ncol = length(defs_names))
  colnames(table.null) <- defs_names

  for (i in 1:R) {
    data_permute <- data
    N_permute <- nrow(data_permute)
    # Sample A1, A2 jointly (without replacement)
    data_A <- data.frame(A1 = data_permute$A1, A2 = data_permute$A2)
    inds_A <- sample(N_permute)
    data_A_sampled <- data_A[inds_A,]

    data_permute$A1 <- data_A_sampled[,1]
    data_permute$A2 <- data_A_sampled[,2]

    # Estimate nuisance parameters and Y0
    est_list_permute <- get_est_analysis(data_permute, cutoff = cutoff, pi_model = NULL,
                                         pi_model_type = pi_model_type, pi_model_seed = pi_model_seed, pi_xvars = pi_xvars)

    ## Get small subgroup unfairness metrics
    if(estimator_type %in% c("small_internal", "small_borrow")) {
      ### Get outcome and P(A=a) models
      models_temp <- get_models_small(data=data_permute, cutoff=cutoff, estimator_type=estimator_type,
                                      outcome_xvars=outcome_xvars, outcome_model_type=outcome_model_type,
                                      pa_xvars_int=pa_xvars_int, fit_method_int=fit_method_int, nfolds=nfolds,
                                      data_external=data_external, fit_method_ext=fit_method_ext, borrow_metric=borrow_metric, pa_xvars_ext=pa_xvars_ext,
                                      pa_model_ext=pa_model_ext)
      ### Get metrics
      defs_permute <- get_defs_analysis(dplyr::mutate(data_permute, Y0 = est_list_permute[['Y0est']], pi = est_list_permute[['pi']]),
                                        cutoff=cutoff, estimator_type=estimator_type,
                                        preds_out = models_temp$preds_out, preds_pa = models_temp$preds_pa)[['defs']]
    } else {
      ## Get original unfairness metrics
      defs_permute <- get_defs_analysis(dplyr::mutate(data_permute, Y0 = est_list_permute[['Y0est']], pi = est_list_permute[['pi']]),
                                        cutoff=cutoff, estimator_type=estimator_type)[['defs']]
    }
    # Add a row to the null distribution table
    table.null[i,] <- defs_permute
  }
  return(table.null)
}

#' Main function: Estimation of nuisance parameters and unfairness metrics.
#'
#' @param data Dataframe or tibble. Must include at least A1, A2, Y (binary 0/1),
#'  D (binary 0/1), covariates used to train propensity score model,
#'  S_prob (probability).
#' @param cutoff Classification cutoff. Must be a single numeric value.
#' @param estimator_type Type of estimation: standard ('standard'), small subgroup internal-only ('small_internal'),
#'  or small subgroup with data borrowing ('small_borrow').
#' @param gen_null T/F: generate null distributions.
#' @param R_null Number of replications for permutation null distribution (default 500).
#' @param bootstrap Obtain bootstrap estimates using rescaled method ('rescaled') or
#'  no bootstrap ('none', the default).
#' @param B Number of bootstrap resamples (default 500).
#' @param m_factor Fractional power for calculating resample size (default 0.75).
#' @param pi_model Pre-fit propensity score model.
#' @param pi_model_type Type of propensity score model: logistic regression ('glm') or random forest ('rf')
#' @param pi_model_seed Numeric vector of random seeds for random forest propensity score model. Required if pi_model_type is 'rf'.
#'  Must be length 2 if pi_model is not specified and length 1 if pi_model is specified.
#'  First element is set prior to model fitting; second element (or only element if pi_model specified) is set prior to prediction.
#' @param pi_xvars Character vector of covariates used in propensity score model.
#' @param outcome_model_type Type of outcome model: logistic regression ('glm'), random forest ('rf'),
#'  or neural network ('neural_net')
#' @param outcome_xvars Character vector of covariates used in outcome models.
#' @param fit_method_int Type of internal data model for P(A=a): multinomial ('multinomial') or neural network ('neural_net')
#' @param nfolds Number of folds for cross-fitting (default 5).
#' @param pa_xvars_int Character vector of covariates used in internal P(A=a) model.
#' @param data_external External data set. Must contain at least A1, A2, and a subset of covariates used to train P(A=a) model.
#' @param fit_method_ext Type of external data model for P(A=a): multinomial ('multinomial') or neural network ('neural_net')
#' @param pa_xvars_ext Character vector of covariates used in external P(A=a) model.
#' @param borrow_metric Metric for data borrowing: AUC ('auc') or Brier score ('brier')
#' @param pa_model_ext Pre-fit external P(A=a) model. Can be either multinomial or neural network model produced by the package nnet.
#'  Model type must match fit_method_ext.
#'
#' @returns List with the following components:
#'  * defs: estimated definitions
#'  * estY0: mean estimated Y0 for each protected group
#'  * table_null: null distribution table (if gen_null = T)
#'  * boot_out: bootstrap estimates for metrics and error rates (if bootstrap = 'rescaled')
#'  * est.choice: Input data frame with Y0 estimates added
#'
#' @export

analysis_estimation <- function(data, cutoff, estimator_type = c('standard', 'small_internal', 'small_borrow'),
                                gen_null = c(T,F), R_null=500, bootstrap = c('none', 'rescaled'), B=500, m_factor = 0.75,
                                pi_model=NULL, pi_model_type = c('glm', 'rf'), pi_model_seed = NULL, pi_xvars,
                                outcome_model_type = c(NULL, 'glm','rf','neural_net'), outcome_xvars=NULL,
                                fit_method_int = c(NULL, 'multinomial', 'neural_net'), nfolds=5, pa_xvars_int=NULL,
                                data_external=NULL, fit_method_ext = c(NULL, 'multinomial', 'neural_net'), pa_xvars_ext=NULL, borrow_metric = c(NULL, 'auc','brier'),
                                pa_model_ext = NULL) {

  # Check inputs
  estimator_type <- match.arg(estimator_type)
  pi_model_type <- match.arg(pi_model_type)
  gen_null <- checkmate::assert_logical(gen_null)
  bootstrap <- match.arg(bootstrap)
  outcome_model_type <- match.arg(outcome_model_type)
  fit_method_int <- match.arg(fit_method_int)
  fit_method_ext <- match.arg(fit_method_ext)
  borrow_metric <- match.arg(borrow_metric)
  cutoff <- checkmate::assert_number(cutoff, lower=0, upper=1)
  nfolds <- checkmate::assert_number(nfolds, lower=1, upper=10)
  pi_xvars <- checkmate::assert_character(pi_xvars)
  outcome_xvars <- checkmate::assert_character(outcome_xvars, null.ok = T)
  pa_xvars_int <- checkmate::assert_character(pa_xvars_int, null.ok = T)
  pa_xvars_ext <- checkmate::assert_character(pa_xvars_ext, null.ok = T)
  m_factor <- checkmate::assert_number(m_factor, lower=0, upper=1)

  if(pi_model_type == "rf") {
    if(is.null(pi_model)) {
      pi_model_seed <- checkmate::assert_numeric(len = 2)
    } else {
      pi_model_seed <- checkmate::assert_numeric(len = 1)
    }
  }
  if(estimator_type=="small_borrow" & ((is.null(data_external)&is.null(pa_model_ext)) | is.null(fit_method_ext) | is.null(borrow_metric) | is.null(pa_xvars_ext))) stop("'fit_method_ext', 'borrow_metric', 'pa_xvars_ext', and either 'data_external' or 'pa_model_ext' must be specified if 'estimator_type' is 'small_borrow'.")
  if(estimator_type%in%c("small_internal", "small_borrow") & (is.null(pa_xvars_int) | is.null(fit_method_int) | is.null(nfolds) | is.null(outcome_xvars) | is.null(outcome_model_type))) stop("'pa_xvars_int', 'fit_method_int', 'nfolds', 'outcome_xvars', 'outcome_model_type' must be specified if 'estimator_type' is 'small_borrow' or 'small_internal'.")

  ## Required data columns and types
  if(!(all(c("A1","A2","Y","D","S_prob") %in% colnames(data)) &
       length(unique(data$Y))==2 & length(unique(data$D))==2)) stop("'data' must contain columns: A1, A2, Y (binary), D (binary), S_prob (continuous), and propensity score model covariates.")
  if(!all(pi_xvars %in% colnames(data))) stop("'data' must contain all propensity score model covariates (pi_xvars).")

  if(estimator_type%in%c("small_internal", "small_borrow")) {
    if(!all(outcome_xvars %in% colnames(data))) stop("'data' must contain all outcome model covariates (outcome_xvars).")
    if(!all(pa_xvars_int %in% colnames(data))) stop("'data' must contain all internal P(A=a) model covariates (pa_xvars_int).")
  }
  if(estimator_type=="small_borrow") {
    if(!all(pa_xvars_ext %in% pa_xvars_int)) stop("All covariates in 'pa_xvars_ext' must be present in 'pa_xvars_int'.")
    if(!all(pa_xvars_ext %in% colnames(data_external))) stop("'data_external' must contain all external P(A=a) model covariates (pa_xvars_ext).")
  }

  ## Convert Y, D levels other than 0/1 to 0/1, make numeric
  if(!identical(sort(levels(as.factor(data$Y))), c("0","1"))) {
    Ytemp <- replace(data$Y, as.factor(data$Y)==levels(as.factor(data$Y))[1], 0)
    Ytemp <- replace(Ytemp, Ytemp==levels(as.factor(data$Y))[2], 1)
    data$Y <- Ytemp
  }
  if(!identical(sort(levels(as.factor(data$D))), c("0","1"))) {
    Dtemp <- replace(data$D, as.factor(data$D)==levels(as.factor(data$D))[1], 0)
    Dtemp <- replace(Dtemp, Dtemp==levels(as.factor(data$D))[2], 1)
    data$D <- Dtemp
  }
  ## Drop non-complete cases and give a warning
  data_complete_ind <- stats::complete.cases(data)
  data_complete <- dplyr::filter(data, data_complete_ind)
  if(nrow(data_complete) != nrow(data)) {
    warning("'data' contains missing values, dropping incomplete rows.")
    data <- data_complete
  }
  ## External data: required columns and types
  if(!is.null(data_external)) {
    if(!all(c("A1","A2") %in% colnames(data_external))) stop("'data' must contain columns: A1, A2, and covariates for P(A=a) model.")
    data_ext_complete_ind <- stats::complete.cases(data_external)
    data_ext_complete <- dplyr::filter(data_external, data_ext_complete_ind)
    if(nrow(data_ext_complete) != nrow(data_external)) {
      warning("External data contains missing values, dropping incomplete rows")
      data_external <- data_ext_complete
    }
    data_external$A1A2 <- as.factor(paste0(data_external$A1, data_external$A2))
  }
  ## pa_xvars_ext must be a subset of pa_xvars_int
  if(!all(pa_xvars_ext %in% pa_xvars_int)) stop("All external P(A=a) covariates must be included in internal covariates ('pa_xvars_ext' must be a subset of 'pa_xvars_int').")

  # Make protected characteristics factors if they aren't already
  data$A1 <- as.factor(as.character(data$A1))
  data$A2 <- as.factor(as.character(data$A2))

  # Set metric names and return vector length
  if(nlevels(data$A2)>1) {
    A1_levels <- levels(data$A1)
    A2_levels <- levels(data$A2)
    A1A2_grid <- expand.grid(A2_levels, A1_levels)
  } else {
    A1A2_grid <- data.frame(Var2 = levels(data$A1), Var1 = 1)
  }
  A1A2_names <- paste0(A1A2_grid[,"Var2"], A1A2_grid[,"Var1"])

  defs_names_standard <- c(paste0("fpr_", A1A2_names), paste0("fnr_", A1A2_names), paste0("cfpr_", A1A2_names), paste0("cfnr_", A1A2_names),
                  paste0("cfpr_marg_A1_", levels(data$A1)), paste0("cfpr_marg_A2_", levels(data$A2)), paste0("cfnr_marg_A1_", levels(data$A1)), paste0("cfnr_marg_A2_", levels(data$A2)),
                  "avg_neg", "max_neg", "var_neg",
                  "avg_pos", "max_pos", "var_pos")
  defs_names_small <-  c(defs_names_standard,
                         paste0("cfpr_small_", A1A2_names), paste0("cfnr_small_", A1A2_names),
                         "avg_neg_small", "max_neg_small", "var_neg_small",
                         "avg_pos_small", "max_pos_small", "var_pos_small")

  if(estimator_type %in% c("small_internal", "small_borrow")) {
    defs_names <- defs_names_small
  } else {
    defs_names <- defs_names_standard
  }

  # If gen_null option is selected, generate the null distribution
  if(gen_null) {
    # Get null distribution
    table.null.temp <- analysis_nulldist(data, R = R_null, cutoff = cutoff,
                                         pi_model_type = pi_model_type, pi_model_seed = pi_model_seed, pi_xvars = pi_xvars,
                                         defs_names = defs_names, estimator_type = estimator_type,
                                         outcome_model_type=outcome_model_type, outcome_xvars=outcome_xvars,
                                         fit_method_int=fit_method_int, nfolds=nfolds, pa_xvars_int=pa_xvars_int,
                                         data_external=data_external, fit_method_ext=fit_method_ext, borrow_metric=borrow_metric, pa_xvars_ext=pa_xvars_ext,
                                         pa_model_ext=pa_model_ext)
  }
  # Get bootstrap estimates: resample conditional on observed protected group proportions
  if (bootstrap == 'rescaled') {
    boot.out <- bs_rescaled_analysis(data, B = B, m_factor = m_factor, cutoff = cutoff,
                                     pi_model_type = pi_model_type, pi_model_seed = pi_model_seed, pi_xvars = pi_xvars,
                                     defs_names = defs_names, estimator_type = estimator_type,
                                     outcome_model_type=outcome_model_type, outcome_xvars=outcome_xvars,
                                     fit_method_int=fit_method_int, nfolds=nfolds, pa_xvars_int=pa_xvars_int,
                                     data_external=data_external, fit_method_ext=fit_method_ext, borrow_metric=borrow_metric, pa_xvars_ext=pa_xvars_ext,
                                     pa_model_ext=pa_model_ext)
  }

  # Estimate nuisance parameters and Y^0
  est_choice_list <- get_est_analysis(data, cutoff = cutoff,
                                      pi_model = pi_model, pi_model_type = pi_model_type, pi_xvars = pi_xvars)

  ## Add estimated Y0 and pi to test data for return
  est_choice_temp <- data
  est_choice_temp$Y0est <- est_choice_list[['Y0est']]
  est_choice_temp$pi <- est_choice_list[['pi']]

  ## Get small subgroup unfairness metrics
  if(estimator_type %in% c("small_internal", "small_borrow")) {
    ### Get outcome and P(A=a) models
    models_temp <- get_models_small(data=data, cutoff=cutoff, estimator_type=estimator_type,
                                    outcome_xvars=outcome_xvars, outcome_model_type=outcome_model_type,
                                    pa_xvars_int=pa_xvars_int, fit_method_int=fit_method_int, nfolds=nfolds,
                                    data_external=data_external, fit_method_ext=fit_method_ext, borrow_metric=borrow_metric, pa_xvars_ext=pa_xvars_ext,
                                    pa_model_ext=pa_model_ext)
    ### Get metrics
    defs_temp <- get_defs_analysis(dplyr::mutate(data, Y0 = est_choice_list[['Y0est']], pi = est_choice_list[['pi']]), cutoff = cutoff,
                                   estimator_type = estimator_type,
                                   preds_out = models_temp$preds_out, preds_pa = models_temp$preds_pa)
  } else {
  ## Get original unfairness metrics
    defs_temp <- get_defs_analysis(dplyr::mutate(data, Y0 = est_choice_list[['Y0est']], pi = est_choice_list[['pi']]), cutoff = cutoff,
                                   estimator_type = estimator_type)
  }
  names(defs_temp[['defs']]) <- defs_names

  return_list <- list(defs = defs_temp[['defs']], est_choice = est_choice_temp)
  if(exists("table.null.temp")) { return_list[['table_null']] = table.null.temp }
  if(exists("boot.out")) { return_list[['boot_out']] = boot.out }
  if(estimator_type == "small_borrow") {
    return_list[['alpha']] <- models_temp$alpha
  }
  return(return_list)
}
