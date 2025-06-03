#' Markov Chain Monte Carlo (MCMC) for TVP-VAR-SV models under dynamic shrinkage priors
#'
#' \code{shrinkDTVPVAR} samples from the joint posterior distribution of the parameters of a
#' TVP-VAR-SV model as described in Cadonna et al. (2020) and returns the MCMC draws. The prior on
#' the VAR coefficients is a dynamic shrinkage prior, as described in Knaus and Frühwirth-Schnatter (2023).
#' The model can be written as:
#' \deqn{Y_t = c_t + \Phi_{1,t} Y_{t-1} + \Phi_{2,t} Y_{t-2} + \cdots + \Phi_{p,t} Y_{t-p} + \epsilon_t}
#' where \eqn{\epsilon_t \sim \mathcal{N}_m(0, \Sigma_t)}.
#'
#' The elements of the VAR coefficients \eqn{\Phi_{i,t}} are assumed to follow component-wise random walks.
#'
#' For further details concerning the algorithms and the model please refer to the papers by Cadonna et al. (2020) and
#' Knaus and Frühwirth-Schnatter (2023).
#' @param y matrix or data frame containing the time series data. The rows correspond to the time
#' points and the columns to the variables.
#' @param p positive integer indicating the number of lags in the VAR model. The default value is 1.
#' @param mod_type character string that reads either \code{"triple"}, \code{"double"} or \code{"ridge"}.
#' Determines whether the triple gamma, double gamma or ridge prior are used for \code{theta_sr} and \code{beta_mean}.
#' The default is "double".
#' @param const logical value indicating whether a constant should be included in the model. The default value is \code{TRUE}.
#' @param niter positive integer, indicating the number of MCMC iterations
#' to perform, including the burn-in. Has to be larger than or equal to \code{nburn} + 2. The default value is 5000.
#' @param nburn non-negative integer, indicating the number of iterations discarded
#' as burn-in. Has to be smaller than or equal to \code{niter} - 2. The default value is \code{round(niter/2)}.
#' @param nthin positive integer, indicating the degree of thinning to be performed. Every \code{nthin} draw is kept and returned.
#' The default value is 1, implying that every draw is kept.
#' @param display_progress logical value indicating whether the progress bar and other informative output should be
#' displayed. The default value is \code{TRUE}.
#' @param TVP_params_beta \emph{optional} named list containing hyperparameter values for the TVP prior of the beta_mean matrix.
#' Not all have to be supplied, with those missing being replaced by the default values. Any list elements that are misnamed will be ignored and a warning will be thrown.
#' Can either be a list of depth 1, which results in all equations having the same hyperparameters, or a list of lists of length \code{m}, where each sub-list contains the
#' hyperparameters for the respective equation. If sub-lists are provided, they can be unnamed, in which case the order of the elements is assumed to be the same as the order of the equations in the data.
#' Alternatively, they can be named \code{eq1}, \code{eq2}, ..., \code{eqm}.
#' The meaning of the hyperparameters is the same as in \code{\link[shrinkTVP]{shrinkDTVP}}, as this function is used to sample the TVP coefficients.
#' A key difference lies in how the adaptive MH arguments for theta and beta_mean are provided. In this function, the adaptive MH arguments are provided as vectors of length 4,
#' where the elements control the MH steps in following order: \code{a_xi}, \code{a_tau}, \code{c_xi}, \code{c_tau}.
#' E.g. if \code{adaptive = c(FALSE, TRUE, TRUE, FALSE)}, then the MH step for \code{a_xi} and \code{c_tau} are non-adaptive, while the MH step for \code{a_tau} and \code{c_xi} are adaptive.
#' Most hyperparameters to be provided are the same as in \code{\link{shrinkTVPVAR}}, however, the following additional ones are required for the dynamic shrinkage prior:
#' \itemize{
#' \item \code{iid}:  logical. The default value is \code{TRUE}.
#' \item \code{a_psi}: numeric vector of length \code{m*p} if const is FALSE and \code{m*p + 1} if const is TRUE. The default value is a vector filled with 0.5.
#' \item \code{c_psi}: numeric vector of length \code{m*p} if const is FALSE and \code{m*p + 1} if const is TRUE. The default value is a vector filled with 0.5.
#' \item \code{a_rho}: positive, real number. The default value is 2.
#' \item \code{b_rho}:positive, real number between 0 and 1. The default value is 0.95.
#' \item \code{alpha_rho}: numeric vector of length \code{m*p} if const is FALSE and \code{m*p + 1} if const is TRUE. The default value is a vector filled with 0.5.
#' \item \code{beta_rho}: numeric vector of length \code{m*p} if const is FALSE and \code{m*p + 1} if const is TRUE. The default value is a vector filled with 3.
#' \item \code{tuning_par_rho}: positive, real number. The default value is 1.
#' \item \code{adaptive_rho}: logical. If \code{TRUE}, the MH step for rho is adaptive, otherwise it is not. The default value is \code{TRUE}.
#' \item \code{target_rate_rho}: positive, real number. The default value is 0.44.
#' \item \code{batch_size_rho}: positive integer. The default value is 50.
#' \item \code{max_adapt_rho}: positive, real number. The default value is 0.01.
#'  }
#' @param TVP_params_sigma \emph{optional} named list containing hyperparameter values for the TVP prior of the Sigma matrix.
#' The structure is the same as for \code{TVP_params_beta}. The default values are the same as for \code{TVP_params_beta}.
#' @return A list of class \code{"shrinkDTVPVAR"} containing:
#' \item{\code{beta}}{an \code{mcmc.tvp.var} object with the VAR coefficient draws.}
#' \item{\code{beta_mean}}{an \code{mcmc.var} object with the beta_mean draws.}
#' \item{\code{theta_sr}}{an \code{mcmc.var} object with the theta_sr draws.}
#' \item{\code{xi2}}{an \code{mcmc.var} object with the xi2 draws.}
#' \item{\code{c_xi}}{an \code{mcmc} object with the c_xi draws.}
#' \item{\code{kappa2}}{an \code{mcmc.var} object with the kappa2 draws.}
#' \item{\code{kappa2_B}}{an \code{mcmc} object with the kappa2_B draws.}
#' \item{\code{a_xi}}{an \code{mcmc} object with the a_xi draws.}
#' \item{\code{tau2}}{an \code{mcmc.var} object with the tau2 draws.}
#' \item{\code{c_tau}}{an \code{mcmc} object with the c_tau draws.}
#' \item{\code{lambda2}}{an \code{mcmc.var} object with the lambda2 draws.}
#' \item{\code{lambda2_B}}{an \code{mcmc} object with the lambda2_B draws.}
#' \item{\code{a_tau}}{an \code{mcmc} object with the a_tau draws.}
#' \item{\code{Sigma}}{an \code{mcmc.tvp.var} object with the covariance matrix draws.}
#' \item{\code{psi}}{an \code{mcmc.tvp.var} object with the psi draws.}
#' \item{\code{rho_p}}{an \code{mcmc.var} object with the rho draws.}
#' \item{\code{pred_objs}}{a list with objects needed for prediction methods.}
#' \item{\code{final_lambda}}{an \code{mcmc.var} with the values of lambda at time T of the dynamic shrinkage process. Used for predicting.}
#' \item{\code{final_lambda_SIGMA}}{an array with the values of lambda of the variance-covariance matrix Sigma at time T of the dynamic shrinkage process.
#' Used for predicting.}
#' \item{\code{rho_p_SIGMA}}{an array with the rho_p values of the variance-covariance matrix Sigma. Used for predicting.}
#' \item{\code{beta_consts}}{a list of \code{mcmc.tvp} objects with the intercept draws (if \code{const} is \code{TRUE}).}
#' \item{\code{psi_consts}}{a list of \code{mcmc.tvp} objects with the psi draws (if \code{const} is \code{TRUE}).}
#' \item{\code{data}}{a list with the original data used for estimation.}
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#'
#' @references
#' Cadonna, A., Frühwirth-Schnatter, S., & Knaus, P. (2020). Triple the Gamma—A Unifying Shrinkage Prior for Variance and Variable Selection in Sparse State Space and TVP Models. \emph{Econometrics}, 8(2), 20.
#'
#' Knaus, P., Bitto-Nemling, A., Cadonna, A., & Frühwirth-Schnatter, S. (2021). Shrinkage in the Time-Varying Parameter Model Framework Using the \code{R} Package \code{shrinkTVP}. \emph{Journal of Statistical Software}, 100(13), 1–32.
#'
#' Knaus, P., & Frühwirth-Schnatter, S. (2023). The Dynamic Triple Gamma Prior as a Shrinkage Process Prior for Time-Varying Parameter Models. \emph{arXiv preprint} arXiv:2312.10487.
#'
#' @seealso \code{\link{TV_heatmap}}, \code{\link{density_plotter}}, \code{\link{state_plotter}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkDTVPVAR(data, p = 2)
#'
#' # Visualize the results
#' plot(res)
#' plot(res$theta_sr)
#'
#' # Change prior to triple gamma
#' res2 <- shrinkDTVPVAR(data, p = 2, mod_type = "triple")
#'
#' # Modify the hyperparameter setup
#' hyperparam <- list(learn_a_xi = FALSE, learn_c_xi = FALSE,
#'                   learn_a_tau = FALSE, learn_c_tau = FALSE,
#'                   a_xi = 0.5, c_xi = 0.5, a_tau = 0.5, c_tau = 0.5)
#'
#' res3 <- shrinkDTVPVAR(data, p = 2, mod_type = "triple",
#'                     TVP_params_beta = hyperparam,
#'                     TVP_params_sigma = hyperparam)
#' }
#'
#' @export
shrinkDTVPVAR <- function(y,
                          p = 1,
                          mod_type = "double",
                          const = TRUE,
                          niter = 5000,
                          nburn = round(niter/2),
                          nthin = 1,
                          display_progress = TRUE,
                          TVP_params_beta = list(),
                          TVP_params_sigma = list()){

  # Input checking (similar to shrinkTVPVAR)
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("y should be a matrix or a data frame.")
  }
  if (!is.character(mod_type) ||
      !(mod_type %in% c("ridge", "double", "triple"))) {
    stop("mod_type should be one of 'ridge', 'double', or 'triple'.")
  }
  if (!is.logical(const) || length(const) != 1) {
    stop("const has to be a single logical value")
  }
  if (!is.numeric(niter) || length(niter) != 1 || niter <= 0 ||
      !is.numeric(nburn) || length(nburn) != 1 || nburn < 0 ||
      !is.numeric(nthin) || length(nthin) != 1 || nthin <= 0) {
    stop("niter, nburn and nthin have to be single, positive integers (with nburn non-negative)")
  }
  if ((niter - nburn) < 2){
    stop("niter has to be larger than or equal to nburn + 2")
  }
  if ((niter - nburn)/2 < nthin){
    stop("nthin can not be larger than (niter - nburn)/2")
  }
  if (!is.logical(display_progress) || length(display_progress) != 1) {
    stop("display_progress must be a single logical value.")
  }
  if (!is.list(TVP_params_beta)) {
    stop("TVP_params_beta should be a list.")
  }
  if (!is.list(TVP_params_sigma)) {
    stop("TVP_params_sigma should be a list.")
  }

  # Extract index and column names for plotting later
  index <- zoo::index(y)[(p + 1):nrow(y)]
  colnames <- colnames(y)
  if (is.null(colnames)){
    colnames <- 1:ncol(y)
  }

  # Force y to a matrix if it is a data frame
  if (is.data.frame(y)){
    y <- as.matrix(y)
  }

  # Create synthetic lagged covariate matrix using mlag
  X <- mlag(y, p)[(p+1):nrow(y), ]
  if (const == TRUE){
    X <- cbind(1, X)
  }
  Y <- y[(p+1):nrow(y), ]
  m <- ncol(Y)  # number of equations
  d <- ncol(X)  # number of covariates
  N <- nrow(Y)  # number of time points

  # Set up TVP hyperparameters (merging user-supplied with defaults)
  TVP_params_default <- list(beta = list(), sigma = list())
  params_default <- list(
    e1 = 0.5,
    e2 = 0.001,
    d1 = 0.5,
    d2 = 0.001,
    beta_a_xi = 10,
    beta_a_tau = 10,
    alpha_a_xi = 5,
    alpha_a_tau = 5,
    beta_c_xi = 2,
    beta_c_tau = 2,
    alpha_c_xi = 5,
    alpha_c_tau = 5,
    a_tuning_par_xi = 1,
    a_tuning_par_tau = 1,
    c_tuning_par_xi = 1,
    c_tuning_par_tau = 1,
    learn_a_xi = TRUE,
    learn_a_tau = TRUE,
    a_eq_c_xi = FALSE,
    a_eq_c_tau = FALSE,
    a_xi = 0.1,
    a_tau = 0.1,
    learn_c_xi = TRUE,
    learn_c_tau = TRUE,
    c_xi = 0.1,
    c_tau = 0.1,
    learn_kappa2_B = TRUE,
    learn_lambda2_B = TRUE,
    kappa2_B = 20,
    lambda2_B = 20,
    Bsigma_sv = 1,
    a0_sv = 5,
    b0_sv = 1.5,
    bmu = 0,
    Bmu = 1,
    adaptive = rep(TRUE, 4),
    target_rates = rep(0.44, 4),
    batch_sizes = rep(50, 4),
    max_adapts = rep(0.01, 4),

    iid = FALSE,
    a_psi = rep(0.5, d),
    c_psi = rep(0.5, d),
    a_rho = 2,
    b_rho = 0.95,
    alpha_rho = rep(0.5, d),
    beta_rho = rep(3, d),
    tuning_par_rho = 1,
    adaptive_rho = TRUE,
    target_rate_rho = 0.44,
    batch_size_rho = 50,
    max_adapt_rho = 0.01
  )

  for(i in 0:(m-1)){
    TVP_params_default[["beta"]][[paste0("eq", i)]] <-
      if (depth(TVP_params_beta) == 2) {
        if (is.null(names(TVP_params_beta))) {
          list_merger(params_default, TVP_params_beta[[i + 1]])
        } else {
          list_merger(params_default, TVP_params_beta[[paste0("eq", i + 1)]])
        }
      } else {
        list_merger(params_default, TVP_params_beta)
      }
    TVP_params_default[["sigma"]][[paste0("eq", i)]] <-
      if (depth(TVP_params_sigma) == 2) {
        if (is.null(names(TVP_params_sigma))) {
          list_merger(params_default, TVP_params_sigma[[i + 1]])
        } else {
          list_merger(params_default, TVP_params_sigma[[paste0("eq", i + 1)]])
        }
      } else  {
        list_merger(params_default, TVP_params_sigma)
      }
  }

  # Overwrite rho_curr_sds with tuning_par_rho
  for(i in 0:(m-1)){
    TVP_params_default[["beta"]][[paste0("eq", i)]]$rho_curr_sds <- rep(
      TVP_params_default[["beta"]][[paste0("eq", i)]]$tuning_par_rho, d)
    TVP_params_default[["sigma"]][[paste0("eq", i)]]$rho_curr_sds <- rep(
      TVP_params_default[["sigma"]][[paste0("eq", i)]]$tuning_par_rho, i)
  }

  # Set up starting values (also passing along required dynamic fields)
  starting_vals <- list(A_st = array(1, c(m, m, N)),
                        D_st = array(1, c(m, m, N)),
                        const = const,
                        lags = p,
                        do_dat_G = FALSE)
  for (i in 1:N){
    starting_vals$A_st[,,i][upper.tri(starting_vals$A_st[,,i])] <- 0
    starting_vals$D_st[,,i][upper.tri(starting_vals$D_st[,,i]) |
                              lower.tri(starting_vals$D_st[,,i])] <- 0
  }
  for(i in 0:(m-1)){
    starting_vals[["beta"]][[paste0("eq", i)]] <-
      list(beta_st = matrix(0, d, N),
           beta_mean_st = rep(0, d),
           theta_sr_st = rep(1, d),
           tau2_st = rep(0.1, d),
           xi2_st = rep(0.1, d),
           kappa2_st = rep(0.1, d),
           lambda2_st = rep(0.1, d),
           kappa2_B_st = 20,
           lambda2_B_st = 20,
           a_xi_st = 0.1,
           a_tau_st = 0.1,
           c_xi_st = 0.1,
           c_tau_st = 0.1,
           d2_st = 1,
           e2_st = 1,
           sv_mu_st = -10,
           sv_phi_st = 0.5,
           sv_sigma2_st = 1,
           sv_latent_st = rep(0, N),
           h0_st = 0,
           sigma2_st = rep(1, N),
           tuning_pars = c(1, 1, 1, 1),
           psi_st = matrix(1, d, N),
           lambda_0_st = rep(1, d),
           lambda_p_st = matrix(1, d, N),
           kappa_p_st = matrix(1, d, N),
           rho_p_st = rep(0.5, d),
           rho_batches_st = matrix(0,
            nrow = TVP_params_default[["beta"]][[paste0("eq", i)]]$batch_size_rho,
            ncol = d),
           rho_curr_sds_st = rep(1, d),
           rho_batch_nrs_st = rep(1, d),
           rho_batch_pos_st = rep(0, d),
           shrink_inter_st = TRUE,
           inter_column_st = 1)
    starting_vals[["sigma"]][[paste0("eq", i)]] <-
      list(beta_st = matrix(0, i, N + 1),
           beta_mean_st = rep(0, i),
           theta_sr_st = rep(1, i),
           tau2_st = rep(0.1, i),
           xi2_st = rep(0.1, i),
           kappa2_st = rep(0.1, i),
           lambda2_st = rep(0.1, i),
           kappa2_B_st = 20,
           lambda2_B_st = 20,
           a_xi_st = 0.1,
           a_tau_st = 0.1,
           c_xi_st = 0.1,
           c_tau_st = 0.1,
           d2_st = 1,
           e2_st = 1,
           sv_mu_st = -10,
           sv_phi_st = 0.5,
           sv_sigma2_st = 1,
           sv_latent_st = rep(0, N),
           h0_st = 0,
           sigma2_st = rep(1, N),
           tuning_pars = c(1, 1, 1, 1),
           psi_st = matrix(1, i, N),
           lambda_0_st = rep(1, i),
           lambda_p_st = matrix(1, i, N),
           kappa_p_st = matrix(1, i, N),
           rho_p_st = rep(0.5, i),
           rho_batches_st = matrix(0,
            nrow = TVP_params_default[["sigma"]][[paste0("eq", i)]]$batch_size_rho,
            ncol = i),
           rho_curr_sds_st = rep(1, i),
           rho_batch_nrs_st = rep(1, i),
           rho_batch_pos_st = rep(0, i),
           shrink_inter_st = TRUE,
           inter_column_st = 1)
  }

  runtime <- system.time({
    suppressWarnings({
      res <- do_shrinkDTVPVAR(Y, X, mod_type, niter, nburn, nthin,
                              display_progress, TVP_params_default, starting_vals)
    })
  })

  if (res$success_vals$success == FALSE){
    stop(paste0("The sampler failed at iteration ",
                res$success_vals$fail_iter,
                " while trying to ",
                res$success_vals$fail, " in equation ", res$success_vals$fail_eq, ". ",
                "Try rerunning the model. If the sampler fails again, try changing the prior to be more informative. ",
                "If the problem still persists, please contact the maintainer: ",
                maintainer("shrinkTVPVAR")))
  } else {
    res$success_vals <- NULL
  }

  if (display_progress == TRUE){
    cat("Timing (elapsed): ", file = stderr())
    cat(runtime["elapsed"], file = stderr())
    cat(" seconds.\n", file = stderr())
    cat(round((niter + nburn)/runtime[3]), "iterations per second.\n\n", file = stderr())
    cat("Converting results to coda objects... ", file = stderr())
  }

  nsave <- floor((niter - nburn) / nthin)

  if (const == TRUE){
    consts_ind <- seq(from = 1, to = m * (m*p) + 1, by = (m*p) + 1)
    params_ind <- setdiff(1:(m * (m*p + 1)), consts_ind)
  } else {
    params_ind <- 1:(m * (m*p))
  }

  res$beta <- aperm(array(res$beta[params_ind, , ], dim = c(m, p, m, N, nsave)), c(3, 1, 2, 4, 5))
  res$Sigma <- aperm(array(res$Sigma, dim = c(m, N, m, nsave)), c(1, 3, 2, 4))
  res$psi <- aperm(array(res$psi[params_ind, , ], dim = c(m, p, m, N, nsave)), c(3, 1, 2, 4, 5))

  for (j in seq_along(res)) {
    attr(res[[j]], "param_name") <- names(res)[j]
    attr(res[[j]], "index") <- index
    attr(res[[j]], "m") <- m
    attr(res[[j]], "p") <- p
    attr(res[[j]], "const") <- const
    attr(res[[j]], "colnames") <- colnames
    if (length(dim(res[[j]])) == 3 && all(dim(res[[j]]) == c(m, 1, nsave))){
      res[[j]] <- as.mcmc(t(drop(res[[j]])))
    }
  }

  class(res$beta) <- "mcmc.tvp.var"
  class(res$Sigma) <- "mcmc.tvp.var"
  class(res$psi) <- "mcmc.tvp.var"

  for (j in seq_along(res)) {
    if (length(dim(res[[j]])) == 3 && all(dim(res[[j]]) == dim(res$theta_sr))){
      class(res[[j]]) <- "mcmc.var"
    }
  }

  if (const == TRUE){
    res$beta_consts <- list()
    res$psi_consts <- list()
    for (i in 1:m){
      res$beta_consts[[i]] <- t(res$beta_const[i, , ])
      class(res$beta_consts[[i]]) <- "mcmc.tvp"
      attr(res$beta_consts[[i]], "index") <- index

      res$psi_consts[[i]] <- t(res$psi_const[i, , ])
      class(res$psi_consts[[i]]) <- "mcmc.tvp"
      attr(res$psi_consts[[i]], "index") <- index
    }
    res$beta_const <- NULL
    res$psi_const <- NULL
  }

  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }

  res$data$Y <- Y
  res$data$X <- X
  attr(res, "class") <- c("shrinkDTVPVAR", "shrinkTVPVAR")
  attr(res, "N") <- N
  attr(res, "m") <- m
  attr(res, "p") <- p
  attr(res, "nsave") <- nsave
  attr(res, "colnames") <- colnames
  attr(res, "index") <- index
  attr(res, "const") <- const
  attr(res, "niter") <- niter
  attr(res, "nburn") <- nburn
  attr(res, "nthin") <- nthin
  attr(res, "TVP_params_beta") <- TVP_params_default$beta
  attr(res, "TVP_params_sigma") <- TVP_params_default$sigma

  return(res)
}
