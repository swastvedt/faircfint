## code to prepare `ex_data_external` dataset goes here
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))
prob_trunc <- function(p) pmax(pmin(p, 0.995),0.005)

# Parameters
params_list <- list(
  # covariate means
  mean_X_pa = c(1,-1,0.5,-0.5,1),
  # coefficients for A
  A1_A2_means = matrix(c(0.3,-.05,0,0,0,-0.2,
                         -0.05,0,.05,0,0.1,0,
                         0.3,0,-.1,0.4,-0.05,-0.2,
                         -0.4,-.1,0,-0.3,0,0.05,
                         0,0,0,0,0.1,-0.4), nrow = 5, ncol = 6, byrow = T)
)

data_gen_external <- function(params, N) {
  # Create local variables from elements of params list
  list2env(params, env = environment())

  # P(A=a) model covariates
  Xmat_pa <- matrix(nrow = N, ncol = length(mean_X_pa))
  for (j in 1:N) {
    Xmat_pa[j,] = mvtnorm::rmvnorm(1, mean = mean_X_pa, sigma = diag(0.5, nrow = length(mean_X_pa)))
  }
  Xmat_pa <- cbind(Xmat_pa)
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

  A1 <- as.factor(as.character(A1))
  A2 <- as.factor(as.character(A2))
  data_list <- list("A1"=A1, "A2"=A2, "X_pa" = Xmat_pa)

  return(data_list)
}

ex_data_external <- as.data.frame(data_gen_external(params_list, N=10000))

# Save external data
usethis::use_data(ex_data_external, overwrite = TRUE)
