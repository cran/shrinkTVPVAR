#' Nicer printing of shrinkTVPVAR objects
#'
#' @param x a \code{shrinkTVPVAR} object.
#' @param ... Currently ignored.
#'
#' @return Called for its side effects and returns invisibly.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
print.shrinkTVPVAR <- function(x, ...){
  ind <- attr(x, "index")
  cat(paste0("Object containing a fitted TVP-VAR model with:\n",
             " - ", formatC(attr(x, "m"), width = 7), " equations\n",
             " - ", formatC(attr(x, "p"), width = 7), " lags\n",
             " - ", formatC(nrow(x$data$Y), width = 7), " timepoints used during estimation, running from ", min(ind), " to ", max(ind), "\n",
             " - ", formatC(attr(x, "niter"), width = 7), " MCMC draws\n",
             " - ", formatC(attr(x, "nburn"), width = 7), " burn-in\n",
             " - ", formatC(attr(x, "nthin"), width = 7), " thinning\n"))
  invisible
}

#' Calculate fitted historical values for an estimated TVP-VAR-SV model
#'
#' Calculates the fitted values for an estimated TVP-VAR-SV model.
#'
#' @param object A \code{shrinkTVPVAR} object
#' @param ... Currently ignored.
#'
#' @return An object of class \code{shrinkTVPVAR_fit}
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' fitted <- fitted(res)
#'
#' # Visualize fitted values
#' plot(fitted)
#' }
#' @family prediction functions
#' @export
fitted.shrinkTVPVAR <- function(object, ...) {
  m <- attr(object, "m")
  p <- attr(object, "p")
  N <- attr(object, "N")
  nsave <- attr(object, "nsave")


  betas <- array(c(aperm(object$beta, c(4, 2, 3, 5, 1))), c(N, m*p, nsave, m))

  offset <- 1
  if (attr(object, "const")) offset <- 2

  fitted <- apply(betas, c(3,4), function(mat) {rowSums(mat * object$data$X[,offset:ncol(object$data$X)])})

  if (attr(object, "const")) {
    # Have to reconstruct original beta_consts
    recons <- array(do.call(c, object$beta_consts), dim = unlist(attributes(object)[c("nsave", "N", "m")]))

    fitted <- fitted + aperm(recons, c(2, 1, 3))
  }

  res_object <- list()
  for(i in 1:m) {
    res_object[[attr(object, "colnames")[i]]] <- t(fitted[,,i])
    attr(res_object[[i]], "class") <- "mcmc.tvp"
  }

  attr(res_object, "colnames") <- attr(object, "colnames")
  attr(res_object, "index") <- attr(object, "index")
  attr(res_object, "class") <- "shrinkTVPVAR_fit"

  return(res_object)
}

#' Draw from posterior predictive density of a fitted TVP-VAR-SV model
#'
#' \code{forecast_shrinkTVPVAR} draws from the posterior predictive distribution of a fitted TVP-VAR-SV model resulting from a call to
#' \code{shrinkTVPVAR}.
#'
#' @param mod an object of class \code{shrinkTVPVAR}, containing the fitted model.
#' @param n.ahead a single, positive integer indicating the forecasting horizon, i.e. how many time-points into the future
#' the posterior predictive distribution should be sampled from. Can not be larger than the number of rows in \code{newdata}.
#'
#' @return The value returned is a list object of class \code{shrinkTVPVAR_forc} containing the samples from the
#' posterior predictive density.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' forc <- forecast_shrinkTVPVAR(res, n.ahead = 4)
#'
#' # Visualize forecast
#' plot(forc)
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family prediction functions
#'
#' @export
forecast_shrinkTVPVAR <- function(mod, n.ahead = 1){

  # Check if mod is of class shrinkTVP
  if (!inherits(mod, "shrinkTVPVAR")) {
    stop("mod has to be an object of class shrinkTVPVAR")
  }

  if (int_input_bad(n.ahead) || n.ahead == 0) {
    stop("n.ahead has to be a single, positive integer")
  }

  # Number of saved posterior draws, number of equations, number of lags
  nsamp <- attr(mod, "nsave")
  m <-  attr(mod, "m")
  p <- attr(mod, "p")
  N <- attr(mod, "N")
  const <- attr(mod, "const")

  # Container for predictions
  y_pred <- array(0, c(n.ahead, m, nsamp))

  # Last values in sample
  curr_A <- mod$pred_objs$final_A
  curr_h <- log(mod$pred_objs$final_D)

  curr_beta <- array(mod$beta[,,,dim(mod$beta)[4],], dim = c(m, m, p, nsamp))

  if (const) {
    # Have to reconstruct original beta_consts
    recons <- array(do.call(c, mod$beta_consts), dim = unlist(attributes(mod)[c("nsave", "N", "m")]))
    recons <- aperm(recons, c(3, 2, 1))
    curr_beta_const <- recons[,dim(recons)[2],]
  }

  for (ti in 1:n.ahead) {

    # Predict SIGMA
    # First predict lower unitriangular matrix
    curr_A <- curr_A + rnorm(length(curr_A), 0, abs(mod$pred_objs$theta_sr_SIGMA))

    # Predict idiosyncratic stochvol equations
    for (eq in 1:m) {
      mean_h <- mod$pred_objs$sv_mu[, eq] + mod$pred_objs$sv_phi[, eq] * (curr_h[, eq] - mod$pred_objs$sv_mu[, eq])
      curr_h[, eq] <- rnorm(nsamp, mean_h, sqrt(mod$pred_objs$sv_sigma2[, eq]))
    }
    curr_D_sr <- array(apply(exp(curr_h/2), 1, diag), c(m, m, nsamp))

    # Predict beta
    offset <- 0

    if (const) {
      curr_beta_const <- curr_beta_const + rnorm(length(curr_beta_const), 0, abs(mod$theta_sr[,1,]))
      y_pred[ti,,] <- curr_beta_const
      offset <- 1
    }

    curr_beta <- curr_beta + rnorm(length(curr_beta), 0, abs(mod$theta_sr[, (1 + offset):(m * p + offset), ]))

    for (curr_p in 1:p){
      # Differentiate between case where the previous y is unknown and where it is known
      curr_prev_t <- ti - curr_p

      if (curr_prev_t <= 0) {
        curr_prev_y <- mod$data$Y[N + curr_prev_t, ]
        y_pred[ti,,] <- y_pred[ti,,] + apply(curr_beta[,, curr_p, ], 3, function(mat) mat %*% curr_prev_y)
      } else {
        curr_prev_y <- y_pred[ti - curr_p,,]
        y_pred[ti,,] <- y_pred[ti,,] + sapply(1:nsamp, function(i) curr_beta[,, curr_p, i] %*% curr_prev_y[,i])
      }
    }

    chol_SIGMA <- array(sapply(1:nsamp, function(i) curr_A[,,i] %*% curr_D_sr[,,i]), c(m, m, nsamp))

    y_pred[ti,, ] <- y_pred[ti,, ] + apply(chol_SIGMA, 3, function(mat) mat %*% rnorm(m))
  }

  res_object <- list()
  for(i in 1:m) {
    name <- paste0(attr(mod, "colnames")[i], "_forc")
    res_object[[name]] <- list()
    res_object[[name]][["y_pred"]] <- t(y_pred[,i,])
    res_object[[name]][["y_orig"]] <- mod$data$Y[,i]
    attr(res_object[[i]], "class") <- c("shrinkTVP_forc", "shrinkTVPVAR_forc_univ")
    attr(res_object[[i]], "index") <- attr(mod, "index")
  }

  class(res_object) <- c("shrinkTVPVAR_forc")
  attr(res_object, "index") <- attr(mod, "index")
  attr(res_object, "colnames") <- attr(mod, "colnames")

  return(res_object)
}
