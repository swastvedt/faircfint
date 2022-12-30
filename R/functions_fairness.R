prob_trunc <- function(p) pmax(pmin(p, 0.995),0.005)

# Get estimated unfairness metrics
#
# 'get_defs_analysis' returns a list of unfairness definitions.
#
# Note: Indices and/or order of definitions in output of this function are
#  used in formatting functions.
#
# @param data Dataframe or tibble. Must include at least A1, A2, Y (binary 0/1),
#  D (binary 0/1), Y0 (estimated), pi (estimated).
# @param cutoff Classification cutoff. Must be a single numeric value.
#
# @returns
# A list with the following components:
#
# * defs: vector of unfairness definitions with 12+(4*nlevels(A1)nlevels(A2))+2*nlevels(A1)+2*nlevels(A2) items

get_defs_analysis <- function(data, cutoff) {
  # Check inputs
  if(!"D"%in%colnames(data) | !"pi"%in%colnames(data)) stop("'data' must include columns D and pi")

  # Grid of values for protected characteristic combinations
  A1_levels <- levels(data$A1)
  A2_levels <- levels(data$A2)
  A1A2_grid <- expand.grid(A1_levels, A2_levels)

  for(i in 1:nrow(A1A2_grid)) {
    varname <- paste0("A1A2_", i)
    data[varname] <- dplyr::if_else(data$A1 == A1A2_grid[i,"Var1"] & data$A2 == A1A2_grid[i,"Var2"],
                             1, 0)
  }
  # Other variables for fairness definitions
  Y0 <- data$Y0
  S_prob <- data$S_prob
  Y <- as.numeric(data$Y)
  data$S <- dplyr::if_else(S_prob >= cutoff, 1, 0)
  S <- data$S
  D <- as.numeric(data$D)
  pi <- data$pi

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
  cfpr_new_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*S*(1-Y)*(1-D)/(1-pi))/mean(data[,varname]*(1-Y)*(1-D)/(1-pi))
  })

  cfnr_vec <- sapply(seq(1, nrow(A1A2_grid)), function(i) {
    varname <- paste0("A1A2_", i)
    mean(data[,varname]*(1-S)*Y0)/mean(data[,varname]*Y0)
  })

  ### Remove groups with 0 or 1 from cfpr_new_vec, cfnr_vec, fpr_vec, fnr_vec
  cfpr_new_vec_dropped <- cfpr_new_vec[cfpr_new_vec != 0 & cfpr_new_vec != 1]
  cfnr_vec_dropped <- cfnr_vec[cfnr_vec != 0 & cfnr_vec != 1]
  fpr_vec_dropped <- fpr_vec[fpr_vec != 0 & fpr_vec != 1]
  fnr_vec_dropped <- fnr_vec[fnr_vec != 0 & fnr_vec != 1]

  # Non-intersectional (marginal) cFPRs and cFNRs
  cfpr_new_marg_vec_A1 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var1))), function(i) {
    varvec <- dplyr::if_else(data$A1 == levels(data$A1)[i], 1, 0)
    mean(varvec*S*(1-Y)*(1-D)/(1-pi))/mean(varvec*(1-Y)*(1-D)/(1-pi))
  })

  cfpr_new_marg_vec_A2 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var2))), function(i) {
    varvec <- dplyr::if_else(data$A2 == levels(data$A2)[i], 1, 0)
    mean(varvec*S*(1-Y)*(1-D)/(1-pi))/mean(varvec*(1-Y)*(1-D)/(1-pi))
  })

  cfnr_marg_vec_A1 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var1))), function(i) {
    varvec <- dplyr::if_else(data$A1 == levels(data$A1)[i], 1, 0)
    mean(varvec*(1-S)*Y0)/mean(varvec*Y0)
  })
  cfnr_marg_vec_A2 <- sapply(seq(1, nrow(dplyr::distinct(A1A2_grid, .data$Var2))), function(i) {
    varvec <- dplyr::if_else(data$A2 == levels(data$A2)[i], 1, 0)
    mean(varvec*(1-S)*Y0)/mean(varvec*Y0)
  })

  ## Non-intersectional deltas
  cdeltaps_new_A1 <- utils::combn(cfpr_new_marg_vec_A1, 2, FUN=function(x) x[1]-x[2])
  cdeltaps_new_A2 <- utils::combn(cfpr_new_marg_vec_A2, 2, FUN=function(x) x[1]-x[2])

  cdeltans_A1 <- utils::combn(cfnr_marg_vec_A1, 2, FUN=function(x) x[1]-x[2])
  cdeltans_A2 <- utils::combn(cfnr_marg_vec_A2, 2, FUN=function(x) x[1]-x[2])

  ## Intersectional deltas
  ### Positive
  cdeltaps_new <- tryCatch({
    utils::combn(cfpr_new_vec_dropped, 2, FUN=function(x) x[1]-x[2])
  }, error=function(cond) {NULL})

  ### Negative
  cdeltans <- tryCatch({
    utils::combn(cfnr_vec_dropped, 2, FUN=function(x) x[1]-x[2])
  }, error=function(cond) {NULL})

  ## Intersectional observational deltas
  ### Negative
  odeltans <- tryCatch({
    utils::combn(fnr_vec_dropped, 2, FUN=function(x) x[1]-x[2])
  }, error=function(cond) {NULL})

  ### Positive
  odeltaps <- tryCatch({
    utils::combn(fpr_vec_dropped, 2, FUN=function(x) x[1]-x[2])
  }, error=function(cond) {NULL})

  ############Get negative (\Delta^-) and positive (\Delta^+) versions of definitions separately
  ### Average
  if(!is.null(cdeltans)) {
    cdelta_avg_neg <- mean(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_avg_neg <- NA_real_ }

  if(!is.null(cdeltaps_new)) {
    cdelta_avg_pos_new <- mean(sapply(cdeltaps_new, abs), na.rm = T)
  } else{ cdelta_avg_pos_new <- NA_real_ }

  ## Max
  if(!is.null(cdeltans)) {
    cdelta_max_neg <- max(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_max_neg <- NA_real_ }

  if(!is.null(cdeltaps_new)) {
    cdelta_max_pos_new <- max(sapply(cdeltaps_new, abs), na.rm = T)
  } else{ cdelta_max_pos_new <- NA_real_ }

  ## 90th percentile of intersectional cdeltas
  cdelta_quant_neg <- tryCatch(
    stats::quantile(sapply(cdeltans, abs), probs = 0.9, na.rm = T),
    error=function(cond) {
      print(cdeltans)
      return(NA) }
  )
  cdelta_quant_pos_new <- tryCatch(
    stats::quantile(sapply(cdeltaps_new, abs), probs = 0.9, na.rm = T),
    error=function(cond) {
      print(cdeltaps_new)
      return(NA) }
  )

  ## Variational
  if(!is.null(cdeltans)) {
    cdelta_var_neg <- stats::var(sapply(cdeltans, abs), na.rm = T)
  } else { cdelta_var_neg <- NA_real_ }

  if(!is.null(cdeltaps_new)) {
    cdelta_var_pos_new <- stats::var(sapply(cdeltaps_new, abs), na.rm = T)
  } else { cdelta_var_pos_new <- NA_real_ }

  ## Observational intersectional
  if(!is.null(odeltans)) {
    odelta_int_neg <- mean(sapply(odeltans, abs), na.rm = T)
  } else { odelta_int_neg <- NA_real_ }

  if(!is.null(odeltaps)) {
    odelta_int_pos <- mean(sapply(odeltaps, abs), na.rm = T)
  } else { odelta_int_pos <- NA_real_ }

  ## Counterfactual non-intersectional (average of marginals)
  cdelta_marg_neg <- mean(c(sapply(cdeltans_A1, abs), sapply(cdeltans_A2, abs)), na.rm = T)
  cdelta_marg_pos_new <- mean(c(sapply(cdeltaps_new_A1, abs), sapply(cdeltaps_new_A2, abs)), na.rm = T)

  ## Definitions
  defs <- c(fpr_vec, fnr_vec, cfpr_new_vec, cfnr_vec,
            cfpr_new_marg_vec_A1, cfpr_new_marg_vec_A2, cfnr_marg_vec_A1, cfnr_marg_vec_A2,
            cdelta_avg_neg, cdelta_max_neg, cdelta_quant_neg, cdelta_var_neg, odelta_int_neg, cdelta_marg_neg,
            cdelta_avg_pos_new, cdelta_max_pos_new, cdelta_quant_pos_new, cdelta_var_pos_new, odelta_int_pos, cdelta_marg_pos_new
  )
  return(list(defs = defs))
}

# Rescaled, stratified bootstrap for unfairness metrics
#
# @param data Dataframe or tibble. Must include at least A1, A2, Y (binary 0/1),
#  D (binary 0/1), covariates used to train propensity score model,
#  S_prob (continuous risk probability).
# @inheritParams get_defs_analysis
# @param B Number of bootstrap replications.
# @param m_factor Fractional power for calculating resample size (default 0.75).
# @param pi_model_type Type of propensity score model. Options are 'glm' (GLM).
#  or 'rf' (random forest)
# @param pi_model_seed Random seed for random forest propensity score model.
# @param f_lasso Formula for GLM propensity score model.
# @param xvars Character vector of covariates used in propensity score model.
# @param defs_length_base Base length for output of get_defs_analysis.
#
# @returns A matrix of definitions, one row for each of the B bootstrap replications.

bs_rescaled_analysis <- function(data, cutoff, B, m_factor = 0.75,
                                 pi_model_type, pi_model_seed = NULL, f_lasso = NULL, xvars,
                                 defs_length_base) {
  ## seed must be specified if PS model type is 'rf'
  if(pi_model_type == "rf" & is.null(pi_model_seed)) stop("'pi_model_seed' must be specified if 'pi_model_type' is 'rf'")

  # Length of output of get_defs_analysis function
  defs_length <- defs_length_base+4*(nlevels(data$A1)*nlevels(data$A2)) + 2*nlevels(data$A1) + 2*nlevels(data$A2)

  # Empty matrix for results
  boot.out <- matrix(nrow = B, ncol = defs_length)
  # m size for m of n bootstrap
  n_rescaled <- floor(nrow(data)^m_factor)

  for (b in 1:B) {
    # Stratified sample conditional on proportions of (A1, A2), Y, S (binary)
    data$A1A2 <- paste0(data$A1, data$A2)
    data$S_char <- dplyr::if_else(data$S_prob >= cutoff, "1", "0")

    data_bs <- as.data.frame(splitstackshape::stratified(data, c("A1A2", "Y", "S_char"), n_rescaled/nrow(data), replace = F))
    # Get Y0 estimator
    ests <- get_est_analysis(data_bs, cutoff = cutoff,
                             pi_model_type = pi_model_type, pi_model_seed = pi_model_seed,
                             f_lasso = f_lasso, xvars = xvars)
    # Get unfairness definitions and mean Y0 in each group
    defs_temp <- get_defs_analysis(dplyr::mutate(data_bs, Y0 = ests[['ipw']][['Y0est']], pi = ests[['ipw']][['pi']]), cutoff = cutoff)
    # Add row to matrix
    boot.out[b,] <- defs_temp[['defs']]
  }
  return(boot.out)
}

# Get estimated propensity score and $Y^0$ weighted estimator
#
# @inheritParams bs_rescaled_analysis
# @param pi_model Propensity score model object.
#
# @returns List with the following components:
#  * Y0est: Vector of inverse probability weighted estimates of $Y^0$.
#  * pi: Vector of estimated propensity scores.
#

get_est_analysis <- function(data, cutoff,
                             pi_model = NULL, pi_model_type, pi_model_seed = NULL, f_lasso = NULL,
                             xvars) {
  # Check inputs
  ## pi_model and type must be specified
  if(is.null(pi_model_type)) stop("'pi_model_type' must be specified")
  if(is.null(pi_model) & pi_model_type == "glm" & is.null(f_lasso)) stop("'f_lasso' must be specified if 'pi_model' is null and type is 'glm'")

  # Blank list for storing results
  est <- list("ipw" = "ipw")
  est_choice <- lapply(est, function(x) list("Y0est" = NULL, "pi" = NULL))
  # Classify predictions according to cutoff
  data$S <- dplyr::if_else(data$S_prob >= cutoff, 1, 0)
  # Data formatted as X matrix
  data_x <- dplyr::select(data, dplyr::all_of(xvars))
  matrix_x <- stats::model.matrix(~.-1, data = data_x)

  ## Get propensity score model if pi_model not specified
  if(is.null(pi_model)) {
    if(pi_model_type == "glm") {
      ### re-calibrate glm model with selected coefficients
      pi_model <- stats::glm(f_lasso, data = data, family = "binomial")

    } else if(pi_model_type == "rf") {
      set.seed(pi_model_seed)
      pi_model <- randomForest::randomForest(matrix_x, as.factor(data$D), importance = T,
                               ntree = 2000, mtry = 90, nodesize = 50)
    }
  }

  # Get propensity score predictions
  ### (bound pi to match truncation of P(D=1) in Mishler paper)
  if(pi_model_type == 'glm') {
    pihat <- prob_trunc(stats::predict(pi_model, newdata = data, type = "response"))
  } else if(pi_model_type == 'rf') {
    pihat <- prob_trunc(stats::predict(pi_model, newdata = matrix_x, type = "prob")[,2])
  } else stop("'pi_model_type' must be 'glm' or 'rf'")
  ## IPW estimates of Y0
  ipwhat <- ((1-data$D)*data$Y)/(1-pihat)
  ## Save Y0 and pi estimates
  est_choice[['ipw']][['Y0est']] <- ipwhat
  est_choice[['ipw']][['pi']] <- pihat

  return(est_choice)
}

# Simulate null distributions.
#
# @inheritParams bs_rescaled_analysis
# @param R Number of replications.
#
# @returns List with the following components:
#  * ipw: Matrix of null unfairness metrics, one row for each replication.

analysis_nulldist <- function(data, R, cutoff, xvars,
                              pi_model_type, pi_model_seed = NULL, f_lasso = NULL,
                              defs_length_base) {
  # Check inputs
  ## seed must be specified if PS model type is 'rf'
  if(pi_model_type == "rf" & is.null(pi_model_seed)) stop("'pi_model_seed' must be specified if 'pi_model_type' is 'rf'")

  # Length of output of get_defs_analysis function
  defs_length <- defs_length_base+4*(nlevels(data$A1)*nlevels(data$A2)) + 2*nlevels(data$A1) + 2*nlevels(data$A2)

  # Set est as placeholder
  est <- list("ipw" = "ipw")

  # List of empty matrices for storing results
  ## Matrix number of columns is length of get_defs_analysis output
  table.null <- lapply(est, function(x) {
    temp_null <- matrix(nrow = R, ncol = defs_length)
    return(temp_null)
  })

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
    est_list_permute <- get_est_analysis(data_permute, cutoff = cutoff,
                                         pi_model_type = pi_model_type, pi_model_seed = pi_model_seed,
                                         f_lasso = f_lasso, xvars = xvars)

    ## Get unfairness definitions
    defs_list_permute <- lapply(est, function(x) {
      defs_temp <- get_defs_analysis(dplyr::mutate(data_permute, Y0 = est_list_permute[[x]][['Y0est']], pi = est_list_permute[[x]][['pi']]), cutoff = cutoff)
      return(defs_temp)
    })

    ## Lists of defs, one list item for each element of 'est'
    defs_permute <- lapply(defs_list_permute, function(x) x[['defs']])

    # Add a row to the null distribution table
    table.null <- lapply(seq_along(table.null), function(x) {
      newdf <- table.null[[x]]
      newdf[i,] <- defs_permute[[x]]
      return(newdf)
    })
  }
  # Un-nest the resulting null distribution table if only one estimator requested
  names(table.null) <- est
  if(length(est) == 1) {table.null <- table.null[[1]]}
  return(table.null)
}


#' Main function: Estimation of nuisance parameters and unfairness metrics.
#'
#' @param data Dataframe or tibble. Must include at least A1, A2, Y (binary 0/1),
#'  D (binary 0/1), covariates used to train propensity score model,
#'  S_prob (continuous risk probability).
#' @param cutoff Classification cutoff. Must be a single numeric value.
#' @param pi_model Propensity score model object.
#' @param pi_model_type Type of propensity score model. Options are 'glm' (GLM).
#'  or 'rf' (random forest)
#' @param pi_model_seed Random seed for random forest propensity score model.
#' @param f_lasso Formula for GLM propensity score model.
#' @param xvars Character vector of covariates used in propensity score model.
#' @param gen_null T/F: generate null distributions.
#' @param R_null Number of replications for null distribution.
#' @param bootstrap Obtain bootstrap estimates using rescaled method ('rescaled') or
#'  no bootstrap ('none', the default).
#' @param m_factor Fractional power for calculating resample size (default 0.75).
#'
#' @returns List with the following components:
#'  * defs: estimated definitions
#'  * estY0: mean estimated Y0 for each protected group
#'  * table_null: null distribution table (if gen_null = T)
#'  * boot_out: bootstrap estimates for metrics and error rates (if bootstrap = 'rescaled')
#'  * est.choice: Input data frame with Y0 estimates added
#'
#' @export

analysis_estimation <- function(data, cutoff,
                                pi_model, pi_model_type, pi_model_seed = NULL, f_lasso = NULL,
                                xvars,
                                gen_null = F, R_null = 1000,
                                bootstrap = 'none', m_factor = 0.75) {
  defs_length_base <- 12
  # Check inputs
  ## bootstrap must be one of 'none', 'rescaled'
  if(!(bootstrap %in% c('none', 'rescaled'))) stop("'bootstrap' must be one of 'none', 'rescaled'")
  ## pi_model and type must be specified
  if(is.null(pi_model) | is.null(pi_model_type)) stop("'pi_model' and 'pi_model_type' must be specified")
  ## seed must be specified if PS model type is 'rf'
  if(pi_model_type == "rf" & is.null(pi_model_seed)) stop("'pi_model_seed' must be specified if 'pi_model_type' is 'rf'")

  # Set names to preserve in output lists
  names(cutoff) <- cutoff
  ## placeholder for 'est'
  est <- list("ipw" = "ipw")
  # Make protected characteristics factors if they aren't already
  data$A1 <- as.factor(as.character(data$A1))
  data$A2 <- as.factor(as.character(data$A2))

  # If gen_null option is selected, generate the null distribution
  if(gen_null) {
    # Get null distribution
    table.null.temp <- analysis_nulldist(data, R = R_null, cutoff = cutoff,
                                         pi_model_type = pi_model_type, pi_model_seed = pi_model_seed,
                                         f_lasso = f_lasso, xvars = xvars, defs_length_base = defs_length_base)
  }
  # Get bootstrap estimates: resample conditional on observed protected group proportions
  if (bootstrap == 'rescaled') {
    boot.out <- lapply(est, function(x) {
      bs_rescaled_analysis(data, B = 1000, m_factor = m_factor, xvars = xvars,
                           cutoff = cutoff,
                           pi_model_type = pi_model_type, pi_model_seed = pi_model_seed,
                           f_lasso = f_lasso, defs_length_base = defs_length_base)
    })
  }

  # Estimate nuisance parameters and Y^0
  # Get a named list of the requested Y0 estimators, one sub-list for each cutoff
  est_choice_list <- lapply(cutoff, function(c) get_est_analysis(data, cutoff = c,
                                                                 pi_model = pi_model, pi_model_type = pi_model_type,
                                                                 xvars = xvars))

  ## Return estimated Y0 vector and test data, one sub-df for each cutoff
  if(length(est) > 1) {
    est_choice_temp <- lapply(seq_along(cutoff), function(ind) {
      df <- data
      for(j in 1:length(est)) {
        colname <- paste0("Y0est_", names(est_choice_list[[ind]][j]))
        colname_pi <- paste0("pihat_", names(est_choice_list[[ind]][j]))
        df[colname] <- est_choice_list[[ind]][[j]][['Y0est']]
        df[colname_pi] <- est_choice_list[[ind]][[j]][['pi']]
      }
      return(df)
    })
  } else {
    est_choice_temp <- lapply(seq_along(cutoff), function(ind) {
      df <- data
      df$Y0est <- est_choice_list[[ind]][[1]][['Y0est']]
      df$pi <- est_choice_list[[ind]][[1]][['pi']]
      return(df)
    })
  }
  names(est_choice_temp) <- cutoff

  ## Get unfairness definitions and mean Y0 in each group, one sub-list for each cutoff
  defs_list <- lapply(seq_along(cutoff), function(ind) {
    lapply(est, function(x) {
      defs_temp <- get_defs_analysis(dplyr::mutate(data, Y0 = est_choice_list[[ind]][[x]][['Y0est']], pi = est_choice_list[[ind]][[x]][['pi']]), cutoff = cutoff[ind])
      return(defs_temp)
    })
  })
  names(defs_list) <- cutoff

  ## Vectors of defs, one sublist for each cutoff, with one set of definitions for each element of 'est'
  defs <- lapply(defs_list, function(l) {
    lapply(l, function(x) x[['defs']])
  })

  # Un-nest definitions and est_choice table if only one cutoff was requested
  if(length(cutoff) == 1) {
    defs <- defs[[1]]
    est_choice_temp <- est_choice_temp[[1]]
    ## Un-nest the definitions further if only one 'est' was requested
    if(length(est) == 1) { defs <- defs[[1]] }
  }

  return_list <- list(defs = defs, est_choice = est_choice_temp)
  if(exists("table.null.temp")) { return_list[['table_null']] = table.null.temp }
  if(exists("boot.out")) { return_list[['boot_out']] = boot.out }
  return(return_list)
}
