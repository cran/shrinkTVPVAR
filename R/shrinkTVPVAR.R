
#' Markov Chain Monte Carlo (MCMC) for TVP-VAR-SV models with shrinkage
#'
#' \code{shrinkTVPVAR} samples from the joint posterior distribution of the parameters of a TVP-VAR-SV
#' model with shrinkage as described in  Cadonna et al. (2020) and returns the MCMC draws.
#' The model can be written as:
#' \deqn{Y_t = c_t + \Phi_{1,t} Y_{t-1} + \Phi_{2,t} Y_{t-2} + \cdots + \Phi_{p,t} Y_{t-p} + \epsilon_t}
#' where \eqn{\epsilon_t \sim \mathcal{N}_m(0, \Sigma_t)}.
#'
#' The elements of the VAR coefficients \eqn{\Phi_{i,t}} are assumed to follow component-wise
#' random walk.
#'
#' For further details concerning the algorithms and the model please refer to the paper by Cadonna et al. (2020).
#'
#' @param y matrix or data frame containing the time series data. The rows correspond to the time
#' points and the columns to the variables.
#' @param p positive integer indicating the number of lags in the VAR model. The default value is 1.
#' @param mod_type character string that reads either \code{"triple"}, \code{"double"} or \code{"ridge"}.
#' Determines whether the triple gamma, double gamma or ridge prior are used for \code{theta_sr} and \code{beta_mean}.
#' The default is "double".
#' @param const logical value indicating whether a constant should be included in the model. The default value is \code{TRUE}.
#' @param niter positive integer, indicating the number of MCMC iterations
#' to perform, including the burn-in. Has to be larger than or equal to \code{nburn} + 2. The default value is 10000.
#' @param nburn non-negative integer, indicating the number of iterations discarded
#' as burn-in. Has to be smaller than or equal to \code{niter} - 2. The default value is \code{round(niter / 2)}.
#' @param nthin positive integer, indicating the degree of thinning to be performed. Every \code{nthin} draw is kept and returned.
#' The default value is 1, implying that every draw is kept.
#' @param display_progress logical value indicating whether the progress bar and other informative output should be
#' displayed. The default value is \code{TRUE}.
#' @param TVP_params_beta \emph{optional} named list containing hyperparameter values for the TVP prior of the beta_mean matrix.
#' Not all have to be supplied, with those missing being replaced by the default values. Any list elements that are misnamed will be ignored and a warning will be thrown.
#' Can either be a list of depth 1, which results in all equations having the same hyperparameters, or a list of lists of length \code{m}, where each sub-list contains the
#' hyperparameters for the respective equation. If sub-lists are provided, they can be unnamed, in which case the order of the elements is assumed to be the same as the order of the equations in the data.
#' Alternatively, they can be named \code{eq1}, \code{eq2}, ..., \code{eqm}.
#' The meaning of the hyperparameters is the same as in \code{\link[shrinkTVP]{shrinkTVP}}, as this function is used to sample the TVP coefficients.
#' A key difference lies in how the adaptive MH arguments are provided. In this function, the adaptive MH arguments are provided as vectors of length 4,
#' where the elements control the MH steps in following order: \code{a_xi}, \code{a_tau}, \code{c_xi}, \code{c_tau}.
#' E.g. if \code{adaptive = c(FALSE, TRUE, TRUE, FALSE)}, then the MH step for \code{a_xi} and \code{c_tau} are non-adaptive, while the MH step for \code{a_tau} and \code{c_xi} are adaptive.
#' The following elements can be supplied:
#' \itemize{
#' \item \code{e1}: positive, real number. The default value is 0.5.
#' \item \code{e2}: positive, real number. The default value is 0.001.
#' \item \code{d1}: positive, real number. The default value is 0.5.
#' \item \code{d2}: positive, real number. The default value is 0.001.
#' \item \code{beta_a_xi}: positive, real number. The default value is 10.
#' \item \code{beta_a_tau}: positive, real number. The default value is 10.
#' \item \code{alpha_a_xi}: positive, real number. The default value is 5.
#' \item \code{alpha_a_tau}: positive, real number. The default value is 5.
#' \item \code{beta_c_xi}: positive, real number. The default value is 2.
#' \item \code{beta_c_tau}: positive, real number. The default value is 2.
#' \item \code{alpha_c_xi}: positive, real number. The default value is 5.
#' \item \code{alpha_c_tau}: positive, real number. The default value is 5.
#' \item \code{a_tuning_par_xi}: positive, real number. The default value is 1.
#' \item \code{a_tuning_par_tau}: positive, real number. The default value is 1.
#' \item \code{c_tuning_par_xi}: positive, real number. The default value is 1.
#' \item \code{c_tuning_par_tau}: positive, real number. The default value is 1.
#' \item \code{learn_a_xi}: logical. The default value is \code{TRUE}.
#' \item \code{learn_a_tau}: logical. The default value is \code{TRUE}.
#' \item \code{a_eq_c_xi}: logical. The default value is \code{FALSE}.
#' \item \code{a_eq_c_tau}: logical. The default value is \code{FALSE}.
#' \item \code{a_xi}: positive, real number. The default value is 0.1.
#' \item \code{a_tau}: positive, real number. The default value is 0.1.
#' \item \code{learn_c_xi}: logical. The default value is \code{TRUE}.
#' \item \code{learn_c_tau}: logical. The default value is \code{TRUE}.
#' \item \code{c_xi}: positive, real number. The default value is 0.1.
#' \item \code{c_tau}: positive, real number. The default value is 0.1.
#' \item \code{learn_kappa2_B}: logical. The default value is \code{TRUE}.
#' \item \code{learn_lambda2_B}: logical. The default value is \code{TRUE}.
#' \item \code{kappa2_B}: positive, real number. The default value is 20.
#' \item \code{lambda2_B}: positive, real number. The default value is 20.
#' \item \code{Bsigma_sv}: positive, real number. The default value is 1.
#' \item \code{a0_sv}: positive, real number. The default value is 5.
#' \item \code{b0_sv}: positive, real number. The default value is 1.5.
#' \item \code{bmu}: real number. The default value is 0.
#' \item \code{Bmu}: positive, real number. The default value is 1.
#' \item \code{adaptive}: logical vector of length 4. The default value is \code{rep(TRUE, 4)}.
#' \item \code{target_rates}: numeric vector of length 4. The default value is \code{rep(0.44, 4)}.
#' \item \code{batch_sizes}: numeric vector of length 4. The default value is \code{rep(50, 4)}.
#' \item \code{max_adapts}: numeric vector of length 4. The default value is \code{rep(0.01, 4)}.
#' }
#' @param TVP_params_sigma \emph{optional} named list containing hyperparameter values for the TVP prior of the Sigma matrix.
#' The structure is the same as for \code{TVP_params_beta}. The default values are the same as for \code{TVP_params_beta}.
#'
#' @return The value returned is a \code{shrinkTVPVAR} object containing:
#' \item{\code{beta}}{an object of class \code{"mcmc.tvp.var"} containing the MCMC draws of the VAR coefficients.}
#' \item{\code{beta_mean}}{an object of class \code{"mcmc.var"} containing the MCMC draws of beta_mean for the VAR coefficients.}
#' \item{\code{theta_sr}}{an object of class \code{"mcmc.var"} containing the MCMC draws of theta_sr for the VAR coefficients.}
#' \item{\code{xi2}}{an object of class \code{"mcmc.var"} containing the MCMC draws of xi2 for the VAR coefficients.}
#' \item{\code{c_xi}}{an object of class \code{"mcmc"} containing the MCMC draws of c_xi for the VAR coefficients.}
#' \item{\code{kappa2}}{an object of class \code{"mcmc.var"} containing the MCMC draws of kappa2 for the VAR coefficients.}
#' \item{\code{kappa2_B}}{an object of class \code{"mcmc"} containing the MCMC draws of kappa2_B for the VAR coefficients.}
#' \item{\code{a_xi}}{an object of class \code{"mcmc"} containing the MCMC draws of a_xi for the VAR coefficients.}
#' \item{\code{tau2}}{an object of class \code{"mcmc.var"} containing the MCMC draws of tau2 for the VAR coefficients.}
#' \item{\code{c_tau}}{an object of class \code{"mcmc"} containing the MCMC draws of c_tau for the VAR coefficients.}
#' \item{\code{lambda2}}{an object of class \code{"mcmc.var"} containing the MCMC draws of lambda2 for the VAR coefficients.}
#' \item{\code{lambda2_B}}{an object of class \code{"mcmc"} containing the MCMC draws of lambda2_B for the VAR coefficients.}
#' \item{\code{a_tau}}{an object of class \code{"mcmc"} containing the MCMC draws of a_tau for the VAR coefficients.}
#' \item{\code{Sigma}}{an object of class \code{"mcmc.tvp.var"} containing the MCMC draws of the covariance matrices.}
#' \item{\code{pred_objs}}{a list containing objects required for prediction methods to work.}
#' \item{\code{beta_consts}}{a list of \code{mcmc.tvp} objects containing the MCMC draws of the intercepts.}
#' \item{\code{data}}{a list containing the input data, as well as the synthetic "covariates" used for estimation.}
#' Note that only the values pertaining to the VAR coefficients are returned.
#' The values for the variance-covariance matrix are not returned.
#'
#' To display the output, use \code{plot} on various elements of the list, as well as the \code{\link{TV_heatmap}},
#' \code{\link{density_plotter}} and \code{\link{state_plotter}} function. Many functions that can be applied to \code{coda::mcmc} objects
#' (e.g. \code{coda::acfplot}) can be applied to all output elements that are \code{coda} compatible.
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#'
#' @seealso \code{\link{TV_heatmap}} \code{\link{density_plotter}} \code{\link{state_plotter}}
#' @references Cadonna, A., Frühwirth-Schnatter, S., & Knaus, P. (2020). "Triple the Gamma—A Unifying Shrinkage Prior for Variance and Variable Selection in Sparse State Space and TVP Models."
#' \emph{Econometrics}, 8(2), 20. <doi:10.3390/econometrics8020020>
#'
#' Knaus, P., Bitto-Nemling, A., Cadonna, A., & Frühwirth-Schnatter, S. (2021) "Shrinkage in the Time-Varying Parameter Model Framework Using the \code{R} Package \code{shrinkTVP}."
#' \emph{Journal of Statistical Software} 100(13), 1–32. <doi:10.18637/jss.v100.i13>
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#'
#' # Visualize the results
#' plot(res)
#'
#' plot(res$theta_sr)
#'
#' # Change prior to triple gamma
#' res2 <- shrinkTVPVAR(data, p = 2, mod_type = "triple")
#'
#' # Modify the hyperparameter setup
#' # Estimate under hierarchical horseshoe prior
#' hyperparam <- list(learn_a_xi = FALSE, learn_c_xi = FALSE,
#'                   learn_a_tau = FALSE, learn_c_tau = FALSE,
#'                   a_xi = 0.5, c_xi = 0.5, a_tau = 0.5, c_tau = 0.5)
#'
#' res3 <- shrinkTVPVAR(data, p = 2, mod_type = "triple",
#'                     TVP_params_beta = hyperparam,
#'                     TVP_params_sigma = hyperparam)
#'
#' # Can also specify different hyperparameters for each equation
#' # gen_TVP_params() is a helper function and returns a
#' # list of lists (if for_each_eq = TRUE) where each sub-list
#' # contains the hyperparameters for the respective equation
#' hyperparam2 <- gen_TVP_params(m = 3, for_each_eq = TRUE)
#'
#' hyperparam2[[1]]$learn_a_xi <- FALSE
#' hyperparam2[[1]]$a_xi <- 0.5
#'
#' res4 <- shrinkTVPVAR(data, p = 2, mod_type = "triple",
#'                      TVP_params_beta = hyperparam2)
#'
#' # Now, a_xi is only fixed in first equation
#' plot(res4$a_xi)
#' }
#'
#' @export
shrinkTVPVAR <- function(y,
                         p = 1,
                         mod_type = "double",
                         const = TRUE,
                         niter = 5000,
                         nburn = round(niter/2),
                         nthin = 1,
                         display_progress = TRUE,
                         TVP_params_beta = list(),
                         TVP_params_sigma = list()){

  # Input checking
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("y should be a matrix or a data frame.")
  }

  if (char_input_bad(mod_type) || !(mod_type %in% c("ridge", "double", "triple"))) {
    stop("mod_type should be one of 'ridge', 'double', or 'triple'.")
  }

  if (bool_input_bad(const)) {
    stop("const has to be a single logical value")
  }

  # Check if all integer inputs are correct
  to_test_int <- c(niter = niter,
                   nburn = nburn,
                   nthin = nthin,
                   p = p)

  bad_int_inp <- sapply(to_test_int, int_input_bad)

  if (any(bad_int_inp)){
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive integer"))

  }

  if ((niter - nburn) < 2){
    stop("niter has to be larger than or equal to nburn + 2")
  }

  if (nthin == 0){
    stop("nthin can not be 0")
  }

  if ((niter - nburn)/2 < nthin){
    stop("nthin can not be larger than (niter - nburn)/2")
  }

  # Check if all boolean inputs are correct
  to_test_bool <- c(const = const,
                    display_progress = display_progress)

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single logical value"))

  }

  if (!is.list(TVP_params_beta)) {
    stop("TVP_params_beta should be a list.")
  }

  if (!is.list(TVP_params_sigma)) {
    stop("TVP_params_sigma should be a list.")
  }



  # Debugging tool
  do_dat_G = FALSE

  # Check inputs
  if (!is.matrix(y) && !is.data.frame(y)) {
    stop("y should be a matrix or a data frame.")
  }

  # Check if all integer inputs are correct
  to_test_int <- c(niter = niter,
                   nburn = nburn,
                   nthin = nthin,
                   p = p)

  bad_int_inp <- sapply(to_test_int, int_input_bad)

  if (any(bad_int_inp)){
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive integer"))

  }

  if (!is.character(mod_type) ||
      !(mod_type %in% c("ridge", "double", "triple"))) {
    stop("mod_type should be one of 'ridge', 'double', or 'triple'.")
  }
  if (!missing(TVP_params_beta) && !is.list(TVP_params_beta)) {
    stop("TVP_params_beta should be a list.")
  }
  if (!missing(TVP_params_sigma) && !is.list(TVP_params_sigma)) {
    stop("TVP_params_sigma should be a list.")
  }

  to_test_bool <- c(const = const,
                    display_progress = display_progress)

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single logical value"))

  }

  # Get index from data, used mainly for plotting methods
  index <- zoo::index(y)[(p + 1):nrow(y)]

  # Get the column names from the data, similarly used for plotting
  colnames <- colnames(y)
  if (is.null(colnames)){
    colnames <- 1:ncol(y)
  }

  # Force to a matrix, as cpp functions expect matrices
  if (is.data.frame(y)){
    y <- as.matrix(y)
  }


  # Use mlag function to correctly lag data to create synthetic "X"
  X <- mlag(y, p)[(p+1):nrow(y), ]

  # Add intercept if user specified it
  if (const == TRUE){
    X <- cbind(1, X)
  }
  Y <- y[(p+1):nrow(y), ]
  # Number of equations
  m <- ncol(Y)
  # Number of "covariates"
  d <- ncol(X)
  # Number of time points used for estimation (length of input - lags)
  N <- nrow(Y)


  # TVP_params contains the hyperparameters, both for the betas and the Sigmas
  TVP_params_default <- list(beta = list(),
                             sigma = list())

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
    max_adapts = rep(0.01, 4))

  for(i in 0:(m-1)){
    TVP_params_default[["beta"]][[paste0("eq", i)]] <-
      # Check if user supplied list of lists or just a list
      if (depth(TVP_params_beta) == 2) {
        # Check if user supplied names for the equations
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

  # starting_vals contains the starting values for the sampled coefficients
  # It is also abused to pass some other parameters to the cpp functions
  starting_vals <- list(A_st = array(1, c(m, m, N)),
                        D_st = array(1, c(m, m, N)),
                        const = const,
                        lags = p,
                        do_dat_G = do_dat_G)

  for (i in 1:N){
    starting_vals$A_st[,,i][upper.tri(starting_vals$A_st[,,i])] <- 0
    starting_vals$D_st[,,i][upper.tri(starting_vals$D_st[,,i]) | lower.tri(starting_vals$D_st[,,i])] <- 0
  }

  for(i in 0:(m-1)){
    starting_vals[["beta"]][[paste0("eq", i)]] <-
      list(beta = matrix(0, d, N + 1),
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
           C0_st = 1,
           sigma2_st = rep(1, N),
           tuning_pars = c(1, 1, 1, 1))

    starting_vals[["sigma"]][[paste0("eq", i)]] <-
      list(beta = matrix(0, i, N + 1),
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
           C0_st = 1,
           sigma2_st = rep(1, N),
           tuning_pars = c(1, 1, 1, 1))


  }



  runtime <- system.time({
    suppressWarnings({
      res <- do_shrinkTVPVAR(Y, X, mod_type, niter, nburn, nthin, display_progress, TVP_params_default, starting_vals)
    })
  })

  # Throw an error if the sampler failed
  if (res$success_vals$success == FALSE){
    stop(paste0("The sampler failed at iteration ",
                res$success_vals$fail_iter,
                " while trying to ",
                res$success_vals$fail, " in equation ", res$success_vals$fail_eq, ". ",
                "Try rerunning the model. ",
                "If the sampler fails again, try changing the prior to be more informative. ",
                "If the problem still persists, please contact the maintainer: ",
                maintainer("shrinkTVPVAR")))
  } else {
    res$success_vals <- NULL
  }

  if (display_progress == TRUE){
    cat("Timing (elapsed): ", file = stderr())
    cat(runtime["elapsed"], file = stderr())
    cat(" seconds.\n", file = stderr())
    cat(round( (niter + nburn) / runtime[3]), "iterations per second.\n\n", file = stderr())
    cat("Converting results to coda objects... ", file = stderr())
  }

  # Number of actually returned samples per parameter
  nsave <- floor((niter - nburn) / nthin)

  # This piece of code finds the indexes of the intercepts and separates them
  # from the other betas.
  # The way the betas are returned from cpp is stacked by equation, i.e.
  # intercept_11, intercept_12, ...
  # betas_11, betas_12, ...
  # intercept_21, intercept_22, ...
  # betas_21, betas_22, ...
  # And so on. This filters out the intercepts (or does not if no intercept is present).
  if (const == TRUE){
    consts_ind <-  seq(from = 1, to = m * (m*p) + 1, by = (m*p) + 1)
    params_ind <- setdiff(1:(m * (m*p + 1)), consts_ind)
  } else {
    params_ind <- 1:(m * (m*p))
  }

  # This brings the betas into the desired output array, with the following order of dimensions:
  # - Phi row
  # - Phi column
  # - lag
  # - time point
  # - sample number
  res$beta <- aperm(array(res$beta[params_ind, , ], dim = c(m, p, m, N, nsave)), c(3, 1, 2, 4, 5))

  # Similarly, this squeezes the samples of Sigma into the following format
  # - Sigma row
  # - Sigma column
  # - time point
  # - sample number
  res$Sigma <- aperm(array(res$Sigma, dim = c(m, N, m, nsave)), c(1, 3, 2, 4))

  # Transform all objects that can be represented as coda objects into coda objects
  # Also pass name of parameter as attr to all objects and some other attributes required for plots to work
  for (j in seq_along(res)) {

    attr(res[[j]], "param_name") <- names(res)[j]
    attr(res[[j]], "index") <- index
    attr(res[[j]], "m") <- m
    attr(res[[j]], "p") <- p
    attr(res[[j]], "const") <- const
    attr(res[[j]], "colnames") <- colnames

    if (length(dim(res[[j]])) == 3 && all(dim(res[[j]]) == c(m, 1, nsave))){
      # This is a vector of parameters, so we need to drop the second dimension
      res[[j]] <- as.mcmc(t(drop(res[[j]])))
    }
  }

  # Transform all objects that can be represented as mcmc.tvp.var objects
  # (i.e. the betas and the Sigmas) into mcmc.tvp.var objects
  class(res$beta) <- "mcmc.tvp.var"
  class(res$Sigma) <- "mcmc.tvp.var"

  # Transform all objects that can be represented as mcmc.var objects
  # (i.e. ) into mcmc.var objects
  for (j in seq_along(res)) {
    if (length(dim(res[[j]])) == 3 && all(dim(res[[j]]) == dim(res$theta_sr))){
      class(res[[j]]) <- "mcmc.var"
    }
  }

  # Lastly, beta_const (if estimated) can be turned into a list of mcmc.tvp objects
  if (const == TRUE){
    res$beta_consts <- list()
    for (i in 1:m){
      res$beta_consts[[i]] <- t(res$beta_const[i, , ])
      class(res$beta_consts[[i]]) <- "mcmc.tvp"
      attr(res$beta_consts[[i]], "index") <- index
    }
    res$beta_const <- NULL
  }

  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }

  # Finally, some attributes and outputs to use in later methods
  res$data$Y <- Y
  res$data$X <- X
  attr(res, "class") <- "shrinkTVPVAR"
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

  return(res)
}
