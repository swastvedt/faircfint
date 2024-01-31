## code to prepare `ex_data_estimation` dataset goes here
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))
prob_trunc <- function(p) pmax(pmin(p, 0.995),0.005)

# Parameters
params_1 <- list(
  nr_marg = 0.5,
  nr_int = 0.4,
  nr_maj = 0.6,
  or_marg = 0.4,
  or_int = 0.6,
  or_maj = 0.2,
  orOpp_maj = 0.3,
  orOpp_margint = 0.2
)
params_list <- list(
  pa = c(0.49, 0.23, 0.13, 0.15),
  p = 4,
  mean_X = c(1, -1, 2, -2),
  z_trt_int = 0.5,
  z_trt_marg = 0.4,
  z_trt_A11_A20 = 0.3,
  z_trt_maj = 0.2,
  # coefficients for P(Y0=1|X)
  alpha_Y0 = logit(params_1$nr_maj),
  beta_Y0 = c(1,1,1,0.75),
  # coefficients for P(Y0=1|A1, A2)
  betaA_Y0 = c(logit(params_1$nr_marg)-logit(params_1$nr_maj), logit(params_1$nr_marg)-logit(params_1$nr_maj), logit(params_1$nr_maj)-2*logit(params_1$nr_marg)+logit(params_1$nr_int)),
  # coefficients for P(D=1|X, Y0=1)
  alpha_D = logit(params_1$or_maj),
  beta_D = c(1,1),
  # coefficients for P(D=1|A1, A2, Y0=1)
  betaA_D = c(logit(params_1$or_marg)-logit(params_1$or_maj), logit(params_1$or_marg)-logit(params_1$or_maj), logit(params_1$or_maj)-2*logit(params_1$or_marg)+logit(params_1$or_int)),
  # coefficients for P(D=1|X, A1, A2, Y0=0)
  alpha_DOpp = logit(params_1$orOpp_maj),
  beta_DOpp = c(1,1),
  betaA_DOpp = c(logit(params_1$orOpp_margint)-logit(params_1$orOpp_maj), logit(params_1$orOpp_margint)-logit(params_1$orOpp_maj), logit(params_1$orOpp_maj)-logit(params_1$orOpp_margint)),
  # Coefficients for adding S to opportunity rate and opposite
  betaS_D = logit(0.1),
  betaS_DOpp = logit(0.1)
)

data_gen <- function(params, type, N, rai = NULL, cutoff) {
  # Create local variables from elements of params list
  list2env(params, env = environment())

  ### A1 and A2 group order:
  #### 1: A1=0, A2=0
  #### 2: A1=1, A2=0
  #### 3: A1=0, A2=1
  #### 4: A1=1, A2=1
  A1_A2 <- t(stats::rmultinom(N, size = 1, prob = c(pa[1], pa[2], pa[3], pa[4])))
  A1_A2_vec <- max.col(A1_A2 != 0, ties.method = 'first')
  A1 <- dplyr::if_else(A1_A2_vec %in% c(2,4), 1, 0)
  A2 <- dplyr::if_else(A1_A2_vec %in% c(3,4), 1, 0)

  X <- matrix(nrow = N, ncol = p)
  for (j in 1:N) {
    X[j,] = mvtnorm::rmvnorm(1, mean = mean_X, sigma = diag(0.3, nrow = p))
  }
  # 1) Probability of having the event under no treatment (need rate)
  Y0 <- stats::rbinom(N, 1,
                 prob = prob_trunc(
                   expit(alpha_Y0 + cbind(X, A1, A2, A1*A2)%*%c(beta_Y0, betaA_Y0))
                 ))
  # 2) Probability of no event under treatment, given event under no treatment (intervention strength)
  Y1 <- stats::rbinom(N, 1, prob = dplyr::case_when(
    A1 == 1 & A2 == 1 & Y0 == 1 ~ 1-z_trt_int,
    A1 == 0 & A2 == 1 & Y0 == 1 ~ 1 - z_trt_marg,
    A1 == 1 & A2 == 0 & Y0 == 1 ~ 1 - z_trt_A11_A20,
    A1 == 0 & A2 == 0 & Y0 == 1 ~ 1 - z_trt_maj,
    TRUE ~ 0
  ))
  # 3) Probability of treatment, by whether or not you would have had the event under no treatment (opportunity rate and its opposite)
  if(type == "pre") {
    or <- alpha_D + cbind(X[,1:2], A1, A2, A1*A2)%*%c(beta_D, betaA_D)
    orOpp <- alpha_DOpp + cbind(X[,1:2], A1, A2, A1*A2)%*%c(beta_DOpp, betaA_DOpp)
  } else if (type == "post") {
    if(!exists("rai")) stop("Specify RAI if type = 'post'")
    new_data <- data.frame("A1" = as.factor(as.character(A1)), "A2" = as.factor(as.character(A2)), "X.1" = X[,1], "X.2" = X[,2], "X.3" = X[,3], "X.4" = X[,4])
    S_prob <- stats::predict(rai, newdata = new_data, type = "prob")[,2]
    S <- dplyr::if_else(S_prob >= cutoff, 1, 0)
    or <- alpha_D + cbind(X[,1:2], A1, A2, A1*A2, S)%*%c(beta_D, betaA_D, betaS_D)
    orOpp <- alpha_DOpp + cbind(X[,1:2], A1, A2, A1*A2, S)%*%c(beta_DOpp, betaA_DOpp, betaS_DOpp)
  } else {stop("'type' must be one of 'pre', 'post'")}

  D <- stats::rbinom(N, 1, prob = prob_trunc(expit(or)))
  Y <- (1-D)*Y0 + D*Y1

  A1 <- as.factor(as.character(A1))
  A2 <- as.factor(as.character(A2))
  data_list <- list("A1"=A1, "A2"=A2, "X"=X, "Y"=Y, "D"=D)

  if(type=="post") {
    data_list[["S"]] <- S
    data_list[["S_prob"]] <- S_prob
  }
  return(data_list)
}

# RAI training data
data_train <- as.data.frame(data_gen(params_list, type = "pre", N = 1000))
data_train$Y <- as.factor(as.character(data_train$Y))

# RAI
rai <- randomForest::randomForest(Y ~ A1 + A2 + X.1 + X.2 + X.3 + X.4, data = data_train, mtry = 6)

# Estimation data
ex_data_estimation <- as.data.frame(data_gen(params_list, type = "post", N=5000, rai = rai, cutoff = 0.5))

# Save estimation data
usethis::use_data(ex_data_estimation, overwrite = TRUE)
