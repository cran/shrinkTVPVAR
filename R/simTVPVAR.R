sim_ar <- function(N, mean = 0, variance, phi = 0.93) {

  c <- mean * (1 - phi)

  cond_var <- variance*(1 - phi^2)

  y <- rep(NA_real_, N)
  y[1] <- rnorm(1, mean, sqrt(variance))

  for (t in 2:N) {
    y[t] <- c + phi*y[t-1] + rnorm(1, 0, sqrt(cond_var))
  }

  return(y)
}

#' Generate synthetic data from a TVP-VAR-SV model
#'
#' \code{simTVPVAR} generates synthetic data from a TVP-VAR-SV model. The data is always generated as to be stationary.
#' This is done via a trial and error approach, where the VAR coefficients are drawn from the data generating process until
#' the VAR process is stationary. As such, very large models might take a long time to generate.
#'
#' @param N integer > 2. Indicates the length of the time series to be
#' generated. The default value is 200.
#' @param p integer > 0. Indicates the number of lags in the VAR model. The default value is 2.
#' @param m integer > 1. Indicates the number of equations in the VAR model. The default value is 3.
#' @param prob_0_beta numeric. Indicates the probability of a zero element in the beta_mean matrix. Can be a single value or a vector of length p.
#' The default value is 0.8.
#' @param prob_0_theta numeric. Indicates the probability of a zero element in the theta matrix. Can be a single value or a vector of length p.
#' The default value is 0.8.
#' @param simsig2_theta_sr numeric. Indicates the standard deviation of the normal distribution from which the elements of the theta matrix are drawn.
#' The default value is 0.2.
#' @param simsig2_beta_mean numeric. Indicates the standard deviation of the normal distribution from which the elements of the beta_mean matrix are drawn.
#' The default value is 0.2.
#' @param intercept logical. Indicates whether an intercept should be included in the model. The default value is TRUE.
#' @param display_progress logical. Indicates whether a progress bar should be displayed. The default value is TRUE.
#'
#' @return The value returned is a list object containing:
#' \itemize{
#'  \item \code{data:} data frame that holds the simulated data.
#'  \item \code{true_vals:} list object containing:
#'   \itemize{
#'    \item \code{Phi:} array containing the true VAR coefficients.
#'    \item \code{Sigma:} array containing the true covariance matrices.
#'    \item \code{theta_sr:} array containing the true standard deviations of the theta matrix.
#'    \item \code{beta_mean:} array containing the true means of the beta matrix.
#'   }
#' }
#'
#' @examples
#' \donttest{
#' # Generate a time series of length 300
#' res <- simTVPVAR(N = 300, m = 3, p = 3)
#'
#' # Estimate a model
#' model <- shrinkTVPVAR(y = res$data, p = 3)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
simTVPVAR <- function(N = 200, p = 2, m = 3, prob_0_beta = 0.8, prob_0_theta = 0.8,
                      simsig2_theta_sr = 0.2, simsig2_beta_mean = 0.2, intercept = TRUE, display_progress = TRUE){

  # Input checking
  to_test_int <- c(N = N,
                   p = p,
                   m = m)

  bad_int_inp <- sapply(to_test_int, int_input_bad)

  if (any(bad_int_inp)){
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive integer"))

  }

  # Check if all numeric inputs are correct
  to_test_num <- list(prob_0_beta = prob_0_beta,
                      prob_0_theta = prob_0_theta,
                      simsig2_theta_sr = simsig2_theta_sr,
                      simsig2_beta_mean = simsig2_beta_mean)


  bad_inp <- sapply(to_test_num, numeric_input_bad)

  if (any(bad_inp)){
    bad_inp_names <- names(to_test_num)[bad_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a real, positive number"))
  }

  if (prob_0_beta < 0 | prob_0_beta > 1){
    stop("prob_0_beta has to be between 0 and 1")
  }

  if (prob_0_theta < 0 | prob_0_theta > 1){
    stop("prob_0_theta has to be between 0 and 1")
  }

  to_test_bool <- c(intercept = intercept,
                    display_progress = display_progress)

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single logical value"))

  }




  # The idea here is that the expected value of the marginal variance is constant
  # across equations
  Sigma <- array(0, c(m, m, N + p))

  mus <- rnorm(m, -5, sqrt(3))

  exp_variance <- exp(mus)

  exp_variance[1] <- exp_variance[1] + (max(exp_variance) - exp_variance[1]) * 1.1

  for (eq in 1:m) {
    Sigma[eq, eq, (p+1):(N + p)] <- exp(stochvol::svsim(N, mu = log(exp_variance)[eq])$latent)
  }


  A <- array(0, c(m, m, N + p))

  left_over_variance <- 0

  for (j in 1:m) {

    left_over_variance <- exp_variance[1] - exp_variance[j]

    for (k in 1:m) {
      if (j > k) {
        if (k == (j - 1)){
          current_chunk <- left_over_variance
        } else {
          current_chunk <- runif(1, 0, left_over_variance)
        }

        curr_variance <- sqrt(current_chunk/exp_variance[k])
        left_over_variance <- left_over_variance - current_chunk

        A[j, k, ] <- sim_ar(N + p, 0, curr_variance, phi = 0.93)
      } else if (j == k) {
        A[j, k, ] = 1
      }
    }
  }

  A_t <- aperm(A, c(2, 1, 3))

  for (t in (p+1):(N+p)) {
    Sigma[,,t] = A[,,t] %*% Sigma[,,t] %*% A_t[,,t]
  }


  theta_srs <- array(NA, dim = c(m, m, p))
  beta_means <- array(NA, dim = c(m, m, p))

  for (i in 1:p){
    currprob <- ifelse(length(prob_0_theta) > 1, rep_len(prob_0_theta, p)[i], prob_0_theta)

    theta_srs[,,i] <- matrix(sample(c(0,1), m*m, T, prob = c(currprob, 1-currprob)), nrow = m, byrow = TRUE)
    theta_srs[,,i] <- theta_srs[,,i] * matrix(rnorm(m*m, 0, sqrt(simsig2_theta_sr)), ncol = m)

    currprob <- ifelse(length(prob_0_beta) > 1, rep_len(prob_0_beta, p)[i], prob_0_beta)

    beta_means[,,i] <- matrix(sample(c(0,1), m*m, T, prob = c(currprob, 1-currprob)), nrow = m, byrow = TRUE)
    beta_means[,,i] <- beta_means[,,i] * matrix(rnorm(m*m, 0, sqrt(simsig2_beta_mean)), ncol = m)
  }

  Phi <- array(NA, dim = c(m, m, p, N + p))

  Phi[,,,1] <- array(rnorm(length(beta_means), beta_means, abs(theta_srs)), dim = c(m, m, p))

  if (display_progress){
    pb <- txtProgressBar(min = 2, max = (N+p), style = 3)
  }
  for (t in 2:(N+p)){

    nonexplosive <- FALSE

    while(nonexplosive == FALSE){
      if (p == 1){
        evol <- matrix(rnorm(length(beta_means), 0, abs(theta_srs)), ncol = m)
      } else {
        evol <- array(rnorm(length(beta_means), 0, abs(theta_srs)), dim = c(m, m, p))
      }
      Phi[,,,t] <- Phi[,,,t-1] + evol

      comp_matrix <- matrix(0, m * p, m * p)
      for (i in 1:p){
        range <- ((i - 1) * m + 1):(i*m)
        comp_matrix[(1:m), range] <- Phi[,,i,t]
      }
      if (p > 1){
        comp_matrix[(m+1):(m*p), 1:(m*(p-1))] <- kronecker(diag(p-1), diag(m))
      }

      eigs <- eigen(comp_matrix)$values

      if (all(sqrt(Re(eigs)^2 + Im(eigs)^2) < 1)){
        nonexplosive <- TRUE
      }

    }

    if (display_progress) setTxtProgressBar(pb, t)
  }

  if (display_progress) close(pb)


  # if (intercept == TRUE){
  #   const <- rnorm(m, 0, 1)
  # } else {
  #   const <- rep(0, m)
  # }



  y <- matrix(0, ncol = N + p, nrow = m)
  y[,1:p] <- matrix(rnorm(m * p, 0, 1), nrow = p)

  for(t in (p+1):(N + p)){

    for (j in 1:p){
      y[,t] <- y[,t] + Phi[,,j,t] %*% y[,t-j]
    }

    y[,t] <- y[,t] + t(chol(Sigma[, , t])) %*% rnorm(m, 0, 1)
  }

  y <- y[,(p+1):(N+p)]
  y <- t(y)

  Phi_out <- array(Phi[,,,(p+1):(N+p)], c(m, m, p, N))
  Sigma_out <- array(Sigma[,,(p+1):(N+p)], c(m, m, N))

  res <- list(data = as.data.frame(y),
              true_vals = list(Phi = Phi_out,
                               Sigma = Sigma_out,
                               theta_sr = theta_srs,
                               beta_mean = beta_means))

  return(res)

}
