data_gen_fixedcf <- function(N) {
  expit <- function(x) 1/(1+exp(-x))

  A1_A2 <- t(stats::rmultinom(N, size = 1, prob = c(0.4, 0.2, 0.2, 0.2)))
  A1_A2_vec <- max.col(A1_A2 != 0, ties.method = 'first')
  A1 <- dplyr::if_else(A1_A2_vec %in% c(2,4), 1, 0)
  A2 <- dplyr::if_else(A1_A2_vec %in% c(3,4), 1, 0)

  X <- matrix(nrow = N, ncol = 4)
  for (j in 1:N) {
    X[j,] = mvtnorm::rmvnorm(1, mean = c(1,-1,2,-2), sigma = diag(0.3, nrow = 4))
  }
  # Probability of having the event under no treatment (need rate)
  Y0 <- stats::rbinom(N, 1,
                      prob = prob_trunc(
                        expit(cbind(X, A1, A2, A1*A2)%*%c(rep(1,4), -1, -0.5, -0.5))
                      ))
  # Probability of no event under treatment, given event under no treatment (intervention strength)
  Y1 <- stats::rbinom(N, 1, prob = dplyr::case_when(
    A1 == 1 & A2 == 1 & Y0 == 1 ~ 1 - 0.3,
    A1 == 0 & A2 == 1 & Y0 == 1 ~ 1 - 0.2,
    A1 == 1 & A2 == 0 & Y0 == 1 ~ 1 - 0.2,
    A1 == 0 & A2 == 0 & Y0 == 1 ~ 1 - 0.2,
    TRUE ~ 0
  ))
  # Risk "prediction"
  S_prob <- rbinom(N, 1, prob = dplyr::case_when(
    Y0 == 1 & A1==0 & A2==0 ~ 0.8,
    Y0 == 1 & A1==1 & A2==0 ~ 0.75,
    Y0 == 1 & A1==0 & A2==1 ~ 0.7,
    Y0 == 1 & A1==1 & A2==1 ~ 0.65,
    Y0 == 0 & A1==0 & A2==0 ~ 0.1,
    Y0 == 0 & A1==1 & A2==0 ~ 0.15,
    Y0 == 0 & A1==0 & A2==1 ~ 0.2,
    Y0 == 0 & A1==1 & A2==1 ~ 0.25,
  ))

  # Opportunity rate and opposite
  or <- 0.1 + cbind(X[,1:2], A1, A2, A1*A2, S_prob)%*%c(1,1, -1,-0.5,-0.25, -0.5)
  orOpp <- 0.2 + cbind(X[,1:2], A1, A2, A1*A2, S_prob)%*%c(1,1, -0.25,-0.5,-1, -0.5)
  D <- stats::rbinom(N, 1, prob = prob_trunc(expit(or)))
  Y <- (1-D)*Y0 + D*Y1

  A1 <- as.factor(as.character(A1))
  A2 <- as.factor(as.character(A2))
  data_list <- list("A1"=A1, "A2"=A2, "X"=X, "Y"=Y, "D"=D, "Y0"=Y0, "S_prob"=S_prob)

  return(data_list)
}
