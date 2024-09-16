#' Generate TVP_params that can be used as input for a TVP-VAR-SV model
#'
#' \code{gen_TVP_params} creates either a list or a list of lists of hyperparameters
#' in the correct format to be used as input for a TVP-VAR-SV model estimated by \code{\link{shrinkTVPVAR}}.
#'
#' @param m The number of equations in the VAR model. Ignored if \code{for_each_eq} is set to FALSE.
#' The default value is 2.
#' @param for_each_eq Logical. If TRUE, a list of lists is returned, where each list contains the
#' hyperparameters for one equation. If FALSE, a single list is returned.
#'
#' @return Either a list containing the hyperparameters for all equations or a list of lists
#' containing the hyperparameters for each equation individually.
#'
#' @examples
#' # For a 5 equation model
#' params <- gen_TVP_params(m = 5)
#'
#' # For a model where all equations share the same hyperparameters
#' params <- gen_TVP_params(for_each_eq = FALSE)
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
gen_TVP_params <- function(m = 2, for_each_eq = TRUE) {

  if (int_input_bad(m) || m == 0) {
    stop("m has to be a single, positive integer")
  }

  if (bool_input_bad(for_each_eq)) {
    stop("for_each_eq has to be a single, logical value")
  }

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

  params_TVP <- list()

  if (for_each_eq) {
    for (i in 1:m) {
      params_TVP[[paste0("eq", i)]] <- params_default
    }

    return(params_TVP)
  } else {
    return(params_default)
  }

  return(params_default)
}

# A function to produce lags
mlag <- function(X, p){
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

# Function that determines depth of a list
depth <- function(this, thisdepth = 0){
  if (!is.list(this) | length(this) == 0){
    return(thisdepth)
  } else {
    return(max(unlist(lapply(this, depth, thisdepth = thisdepth + 1))))
  }
}

# Get some utilities from shrinkTVP
util_list <- c("is.scalar", "numeric_input_bad", "numeric_input_bad_zer",
               "numeric_input_bad_", "int_input_bad", "bool_input_bad",
               "char_input_bad", "lty_input_bad", "list_merger")


for (i in seq_along(util_list)) {
  tmp <- utils::getFromNamespace(util_list[i], "shrinkTVP")
  assign(util_list[i], tmp)
}
