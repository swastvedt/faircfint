## code to prepare `ex_data_small` dataset goes here
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))
prob_trunc <- function(p) pmax(pmin(p, 0.995),0.005)

# Parameters
params_list <- list(
  # covariate means
  mean_X_ps = c(1, -1, 2, -2),
  mean_X_outcome = c(0.5,0.5,-0.5,-0.5,1,-1),
  mean_X_pa = c(1,-1,0.5,-0.5,1),
  # coefficients for A
  A1_A2_means = matrix(c(0.3,-.05,0,0,0,-0.2,
                         -0.05,0,.05,0,0.1,0,
                         0.3,0,-.1,0.4,-0.05,-0.2,
                         -0.4,-.1,0,-0.3,0,0.05,
                         0,0,0,0,0.1,-0.4,
                         0.3,0.1,0.1,-0.05,-0.2,-0.2,
                         0.4,-0.1,-0.1,0.05,-0.05,-0.3), nrow = 7, ncol = 6, byrow = T),
  # P(Y1=a | A1, A2, Y0=1)
  z_trt = c(0.5, 0.4, 0.3, 0.2, 0.2, 0.2),
  # coefficients for P(Y0=1|X)
  alpha_Y0 = logit(0.4),
  beta_Y0 = rep(1,8),
  # coefficients for P(Y0=1|A1, A2)
  betaA_Y0 = logit(c(0.5, 0.4, 0.4, 0.6, 0.4, 0.3)),
  # coefficients for P(D=1|X, Y0=1)
  alpha_D = logit(0.6),
  beta_D = rep(0.7,4),
  # coefficients for P(D=1|A1, A2, Y0=1)
  betaA_D = logit(c(0.4, 0.6, 0.6, 0.2, 0.3, 0.3)),
  # coefficients for P(D=1|X, A1, A2, Y0=0)
  alpha_DOpp = logit(0.3),
  beta_DOpp = rep(0.2,4),
  betaA_DOpp = logit(c(0.3, 0.2, 0.2, 0.2, 0.3, 0.3)),
  # Coefficients for adding S to opportunity rate and opposite
  betaS_D = logit(0.1),
  betaS_DOpp = logit(0.1)
)

data_gen_small <- function(params, type, N, rai = NULL, cutoff) {
  # Create local variables from elements of params list
  list2env(params, env = environment())

  # Propensity score model covariates
  Xmat_ps <- matrix(nrow = N, ncol = length(mean_X_ps))
  for (j in 1:N) {
    Xmat_ps[j,] = mvtnorm::rmvnorm(1, mean = mean_X_ps, sigma = diag(0.3, nrow = length(mean_X_ps)))
  }
  p_ps <- ncol(Xmat_ps)
  # Outcome model covariates
  Xmat_outcome <- matrix(nrow = N, ncol = length(mean_X_outcome))
  for (j in 1:N) {
    Xmat_outcome[j,] = mvtnorm::rmvnorm(1, mean = mean_X_outcome, sigma = diag(0.5, nrow = length(mean_X_outcome)))
  }
  Xmat_outcome <- cbind(Xmat_outcome, Xmat_ps[,2:3])
  p_outcome <- ncol(Xmat_outcome)
  # P(A=a) model covariates
  Xmat_pa <- matrix(nrow = N, ncol = length(mean_X_pa))
  for (j in 1:N) {
    Xmat_pa[j,] = mvtnorm::rmvnorm(1, mean = mean_X_pa, sigma = diag(0.5, nrow = length(mean_X_pa)))
  }
  Xmat_pa <- cbind(Xmat_pa, Xmat_ps[,1], Xmat_outcome[,1])
  p_pa <- ncol(Xmat_pa)
  # Generate A1, A2
  A1_A2_vec <- vector(length = N)
  pa_df <- matrix(nrow = N, ncol = length(A1_A2_means))

  for (j in 1:N) {
    p_vec <- expit(Xmat_pa[j,]%*%A1_A2_means)
    pa_df[j,] <- p_vec
    class_val <- t(rmultinom(1,size=1, prob=p_vec/sum(p_vec)))
    A1_A2_vec[j] <- max.col(class_val != 0, ties.method = 'first')
  }
  A1 <- dplyr::if_else(A1_A2_vec %in% c(2,4,6), 1, 0)
  A2 <- dplyr::case_when(
    A1_A2_vec %in% c(1,2) ~ 0,
    A1_A2_vec %in% c(3,4) ~ 1,
    A1_A2_vec %in% c(5,6) ~ 2)
  Amat <- model.matrix(~A-1, data.frame(A=as.character(A1_A2_vec)))

  # Probability of having the event under no treatment (need rate)
  Y0 <- stats::rbinom(N, 1,
                      prob = prob_trunc(
                        expit(alpha_Y0 + cbind(Xmat_outcome, Amat)%*%c(beta_Y0, betaA_Y0))
                      ))
  # Probability of no event under treatment, given event under no treatment (intervention strength)
  Y1 <- stats::rbinom(N, 1, prob = dplyr::case_when(
    A1 == 1 & A2 == 2 & Y0 == 1 ~ 1- z_trt[6],
    A1 == 0 & A2 == 2 & Y0 == 1 ~ 1 - z_trt[5],
    A1 == 1 & A2 == 1 & Y0 == 1 ~ 1- z_trt[4],
    A1 == 0 & A2 == 1 & Y0 == 1 ~ 1 - z_trt[3],
    A1 == 1 & A2 == 0 & Y0 == 1 ~ 1 - z_trt[2],
    A1 == 0 & A2 == 0 & Y0 == 1 ~ 1 - z_trt[1],
    TRUE ~ 0
  ))
  # Probability of treatment, by whether or not you would have had the event under no treatment (opportunity rate and its opposite)
  if(type == "pre") {
    or <- alpha_D + cbind(Xmat_ps, Amat)%*%c(beta_D, betaA_D)
    orOpp <- alpha_DOpp + cbind(Xmat_ps, Amat)%*%c(beta_DOpp, betaA_DOpp)
  } else if (type == "post") {
    if(!exists("rai")) stop("Specify RAI if type = 'post'")
    new_data <- data.frame("A1" = as.factor(as.character(A1)), "A2" = as.factor(as.character(A2)), "X_outcome.1" = Xmat_outcome[,1], "X_outcome.2" = Xmat_outcome[,2], "X_outcome.3" = Xmat_outcome[,3], "X_outcome.4" = Xmat_outcome[,4],
                           "X_outcome.5" = Xmat_outcome[,5], "X_outcome.6" = Xmat_outcome[,6], "X_outcome.7" = Xmat_outcome[,7], "X_outcome.8" = Xmat_outcome[,8])
    S_prob <- stats::predict(rai, newdata = new_data, type = "prob")[,2]
    S <- dplyr::if_else(S_prob >= cutoff, 1, 0)
    or <- alpha_D + cbind(Xmat_ps, Amat, S)%*%c(beta_D, betaA_D, betaS_D)
    orOpp <- alpha_DOpp + cbind(Xmat_ps, Amat, S)%*%c(beta_DOpp, betaA_DOpp, betaS_DOpp)
  }

  D <- stats::rbinom(N, 1, prob = prob_trunc(expit(or)))
  Y <- (1-D)*Y0 + D*Y1

  A1 <- as.factor(as.character(A1))
  A2 <- as.factor(as.character(A2))
  data_list <- list("A1"=A1, "A2"=A2, "X_pa" = Xmat_pa)

  if(type%in%c("post", "pre")) {
    data_list[["X"]] <- Xmat_ps
    data_list[["X_outcome"]] <- Xmat_outcome
    data_list[["X_pa"]] <- Xmat_pa
    data_list[["Y"]] <- Y
    data_list[["D"]] <- D
  }
  if(type=="post") {
    data_list[["S"]] <- S
    data_list[["S_prob"]] <- S_prob
  }
  return(data_list)
}

# RAI training data
data_train <- as.data.frame(data_gen_small(params_list, type = "pre", N = 1000))
data_train$Y <- as.factor(as.character(data_train$Y))

# RAI
rai <- randomForest::randomForest(Y ~ A1 + A2 + X_outcome.1 + X_outcome.2 + X_outcome.3 + X_outcome.4 + X_outcome.5 + X_outcome.6 + X_outcome.7 + X_outcome.8, data = data_train, mtry = 6)

# Estimation data
ex_data_small <- as.data.frame(data_gen_small(params_list, type = "post", N=5000, rai = rai, cutoff = 0.5))


usethis::use_data(ex_data_small, overwrite = TRUE)
