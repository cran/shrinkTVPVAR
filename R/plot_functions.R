plot_mcmc_tmp <- utils::getFromNamespace("plot.mcmc.tvp", "shrinkTVP")

#' Graphical summary of posterior distribution for a time-varying coefficient matrix in a TVP-VAR model
#'
#' \code{plot.mcmc.tvp} plots empirical posterior quantiles for a time-varying parameter coefficient matrix in a TVP-VAR model.
#'
#' @param x \code{mcmc.tvp.var} object
#' @param lag single integer value, indicating the lag of the time-varying VAR to be plotted. The default value is 1.
#' @param mgp vector of length 3, determining the margin line (in \code{\link[graphics]{par}}) for the plot. The default value is \code{c(1.5, 0.5, 0)}.
#' See \code{\link[graphics]{par}} for more information.
#' @param ylim numeric vector of length 2, determining the y-axis limits of the plot.
#' If missing, the limits are determined by the lowest and largest quantiles of the data.
#' @param ylabs character vector of length m, determining the y-axis labels of the plot.
#' If missing, the labels are taken from the column names of the data.
#' @param mains character vector of length m, determining the main titles of the plot.
#' If missing, the titles are taken from the column names of the data.
#' @param h_borders numeric vector of length 2, determining the horizontal borders of the plot.
#' The first value is the space between the plot and the left border,
#' the second value is the space between the plot and the right border.
#' Both are fractions of the total width of the plot. The default value is \code{c(0.075, 0.05)}.
#' @param w_borders numeric vector of length 2, determining the vertical borders of the plot.
#' The first value is the space between the plot and the top border,
#' the second value is the space between the plot and the bottom border.
#' Both are fractions of the total height of the plot. The default value is \code{c(0.05, 0.05)}.
#' @param ... further arguments to be passed to \code{\link[shrinkTVP]{plot.mcmc.tvp}} (see shrinkTVP package).
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' plot(res$beta)
#'
#' # Plot second lag
#' plot(res$beta, lag = 2)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
state_plotter <- function(x, lag = 1, mgp = c(1.5,0.5,0), ylim, ylabs, mains,
                          h_borders = c(0.075, 0.05), w_borders = c(0.05, 0.05), ...){

  # Get some attributes
  m <- attr(x, "m")
  param_name <- attr(x, "param_name")
  colnames <- attr(x, "colnames")
  p <- attr(x, "p")

  args <- list(...)

  # Input checks
  if (!inherits(x, "mcmc.tvp.var")) {
    stop("x must be of class mcmc.tvp.var")
  }

  if (int_input_bad(lag) | lag > p){
    stop(paste0("lag must be a positive integer, not larger than lag length of model (i.e. ", p, ")"))
  }


  if ((length(mgp) != 3) | any(sapply(mgp, numeric_input_bad_zer))) {
    stop("mgp must be a vector of length 3 with non-negative numeric values")
  }

  if (!missing(ylim)){
    if ((length(ylim) != 2) | any(sapply(ylim, numeric_input_bad))) {
      stop("ylim must be a vector of length 2 with numeric values")
    }
  }

  if (!missing(ylabs)){
    if ((length(ylabs) != m) | any(sapply(ylabs, char_input_bad))) {
      stop("ylabs must be a vector of length m with character values")
    }
  }

  if (!missing(mains)){
    if ((length(mains) != m) | any(sapply(mains, char_input_bad))) {
      stop(paste0("mains must be a vector equal to number of length ", m, " with character values"))
    }
  }

  if (any(sapply(h_borders, numeric_input_bad_zer)) | sum(h_borders) >= 1 | length(h_borders) > 2 |
      (length(h_borders) == 1 & h_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (any(sapply(w_borders, numeric_input_bad_zer)) | sum(w_borders) >= 1 | length(w_borders) > 2 |
      (length(w_borders) == 1 & w_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (length(h_borders) == 1){
    h_borders <- rep(h_borders, 2)
  }

  if (length(w_borders) == 1){
    w_borders <- rep(w_borders, 2)
  }



  # Use custom quantiles if provided, otherwise use default
  if ("probs" %in% names(args)){
    probs <- args$probs
  } else {
    probs <- c(0.025, 0.25, 0.75, 0.975)
  }

  # If ylim is missing, calculate it based on quantiles
  if (missing (ylim)){

    # Have to differentiate between beta and Sigma, as structure is different
    if ((param_name == "beta")){
      quants <- apply(x[,,lag,,], c(1, 2, 3), quantile, c(max(probs), min(probs)))
      Sigma <- FALSE
    } else {
      quants <- apply(x, c(1, 2, 3), quantile, c(max(probs), min(probs)))
      Sigma <- TRUE
    }

    max_y <- max(quants[1,,,])
    min_y <- min(quants[2,,,])

    ylim <- c(min_y, max_y)
  }


  # If ylabs or mains are missing, use colnames
  if (missing(ylabs)) {
    ylabs <- colnames
  } else {
    ylabs <- rep(ylabs, length.out = m)
  }

  if (missing(mains)) {
    mains <- colnames
  } else {
    mains <- rep(mains, length.out = m)
  }


  # Set up the layout
  # Save previous user defined par attributes
  prev_par <- par(no.readonly = TRUE)
  on.exit(par(prev_par))

  par(mar = rep(0, 4), mgp = mgp)

  # Calculate the width and height of each plot
  # Essentially divide up space equally with between 1 - sum(padding)
  width_leftover <- 1 - sum(h_borders)
  curr_widths <- c(h_borders[1], rep(width_leftover/m, m), h_borders[2])

  height_leftover <- 1 - sum(w_borders)
  curr_heights <- c(w_borders[1], rep(height_leftover/m, m), w_borders[2])

  layout(tmp <- matrix(c(rep(1, m + 2), 2:((m+2) * (m) + 1), rep(((m+2) * (m) + 2), m + 2)),
                       nrow = m + 2, ncol = m + 2, byrow = TRUE),
         widths = curr_widths, heights = curr_heights)

  # Loop through the layout and plot the data
  real_ind <- 0L
  for (ind in seq_len(max(tmp))){

    if (ind == 1){
      plot.new()
    } else if (ind %% (m + 2) %in% c(2, 1)){
      plot.new()
    } else if (ind == max(tmp)){
      plot.new()
    } else {
      real_ind <- real_ind + 1L
      # i current row, j current column
      i <- rep(1:m, each = m)[real_ind]
      j <- rep(1:m, m)[real_ind]

      if (Sigma == FALSE | (Sigma == TRUE & j <= i)){
        if (Sigma == FALSE){
          curr_main <- ifelse(i == 1, mains[j], "")
        } else {
          curr_main <- ifelse(i == j, mains[j], "")
        }


        curr_ylab <- ifelse(j == 1, ylabs[i], "")

        curr_xaxt <- ifelse(i == m, "s", "n") # Only last row gets x axis
        curr_yaxt <- ifelse(j == 1, "s", "n") # Only first column gets y axis

        if (Sigma == FALSE){
          dat <- t(x[i, j, lag, , ])
        } else {
          dat <- t(x[i, j, , ])
        }

        attr(dat, "index") <- attr(x, "index")
        plot_mcmc_tmp(dat, ylim = ylim, xaxt = curr_xaxt, yaxt = "n", ...)

        if (i %% 2 == 1 && j == 1){
          axis(2)

        } else if (i %% 2 == 0 && j == m){
          axis(4)
        }

        par(xpd = NA)
        title(main = curr_main, line = 0.5)
        title(ylab = curr_ylab, line = 1.3)
        par(xpd = prev_par$xpd)

      } else {
        plot.new()
      }

    }

  }

  # Reset the plotting parameters
  par(prev_par)

}

#' Kernel density plots of posterior distribution for hyperparameters of time-varying coefficient matrix in a TVP-VAR model
#'
#' \code{density_plotter} plots empirical kernel density estimates of the posterior distribution for hyperparameters
#' of a time-varying parameter coefficient matrix in a TVP-VAR model. \code{beta_mean} and \code{theta_sr} will most
#' likely be of interest here.
#'
#'
#' @param x \code{mcmc.var} object
#' @param lag single integer value, indicating the lag of the time-varying VAR to be plotted. The default value is 1.
#' @param mgp vector of length 3, determining the margin line (in \code{\link[graphics]{par}}) for the plot. The default value is \code{c(1.5, 0.5, 0)}.
#' See \code{\link[graphics]{par}} for more information.
#' @param ylim numeric vector of length 2, determining the y-axis limits of the plot.
#' If missing, the limits are determined by the \code{density} function.
#' @param xlim numeric vector of length 2, determining the x-axis limits of the plot.
#' If missing, the limits are determined by the minimum and maximum values of the data.
#' @param ylabs character vector of length m, determining the y-axis labels of the plot.
#' If missing, the labels are taken from the column names of the data.
#' @param mains character vector of length m, determining the main titles of the plot.
#' If missing, the titles are taken from the column names of the data.
#' @param h_borders numeric vector of length 2, determining the horizontal borders of the plot.
#' The first value is the space between the plot and the left border,
#' the second value is the space between the plot and the right border.
#' Both are fractions of the total width of the plot. The default value is \code{c(0.075, 0.05)}.
#' @param w_borders numeric vector of length 2, determining the vertical borders of the plot.
#' The first value is the space between the plot and the top border,
#' the second value is the space between the plot and the bottom border.
#' Both are fractions of the total height of the plot. The default value is \code{c(0.05, 0.05)}.
#' @param ... further arguments to be passed to \code{\link[graphics]{plot}}.
#'
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' plot(res$theta_sr)
#'
#' # Plot second lag
#' plot(res$theta_sr, lag = 2)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
density_plotter <- function(x, lag = 1, mgp = c(1.5,0.5,0), ylim, xlim, ylabs, mains,
                            h_borders = c(0.075, 0.05), w_borders = c(0.05, 0.05), ...) {
  # Get some attributes
  p <- attr(x, "p")
  const <- attr(x, "const")
  m <- attr(x, "m")
  colnames <- attr(x, "colnames")

  # Input checks
  if (!inherits(x, "mcmc.var")) {
    stop("x must be of class mcmc.var")
  }

  if (int_input_bad(lag) | lag > p){
    stop(paste0("lag must be a positive integer, not larger than", p))
  }

  if ((length(mgp) != 3) | any(sapply(mgp, numeric_input_bad_zer))) {
    stop("mgp must be a vector of length 3 with non-negative numeric values")
  }

  if (!missing(ylim)){
    if ((length(ylim) != 2) | any(sapply(ylim, numeric_input_bad))) {
      stop("ylim must be a vector of length 2 with numeric values")
    }
  }

  if (!missing(ylabs)){
    if ((length(ylabs) != m) | any(sapply(ylabs, char_input_bad))) {
      stop("ylabs must be a vector of length m with character values")
    }
  }

  if (!missing(mains)){
    if ((length(mains) != m) | any(sapply(mains, char_input_bad))) {
      stop(paste0("mains must be a vector of length ", m, " with character values"))
    }
  }

  if (any(sapply(h_borders, numeric_input_bad_zer)) | sum(h_borders) >= 1 | length(h_borders) > 2 |
      (length(h_borders) == 1 & h_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (any(sapply(w_borders, numeric_input_bad_zer)) | sum(w_borders) >= 1 | length(w_borders) > 2 |
      (length(w_borders) == 1 & w_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (length(h_borders) == 1){
    h_borders <- rep(h_borders, 2)
  }

  if (length(w_borders) == 1){
    w_borders <- rep(w_borders, 2)
  }


  if (missing(xlim)) {
    xlim <- c(min(x), max(x))
  }

  if (missing(ylabs)) {
    ylabs <- colnames
  } else {
    ylabs <- rep(ylabs, length.out = m)
  }

  if (missing(mains)) {
    mains <- colnames
  } else {
    mains <- rep(mains, length.out = m)
  }

  prev_par <- par(no.readonly = TRUE)
  on.exit(par(prev_par))
  par(mar = rep(0, 4), mgp = mgp)

  width_leftover <- 1 - sum(h_borders)
  curr_widths <- c(h_borders[1], rep(width_leftover/m, m), h_borders[2])

  height_leftover <- 1 - sum(w_borders)
  curr_heights <- c(w_borders[1], rep(height_leftover/m, m), w_borders[2])

  layout(tmp <- matrix(c(rep(1, m + 2), 2:((m+2) * (m) + 1), rep(((m+2) * (m) + 2), m + 2)), nrow = m + 2, ncol = m + 2, byrow = TRUE),
         widths = curr_widths, heights = curr_heights)

  # Idea of this function very similar to the state plotter, but plotting densities instead
  real_ind <- 0L
  for (ind in seq_len(max(tmp))){

    if (ind == 1){
      plot.new()
    } else if (ind %% (m + 2) %in% c(2, 1)){
      plot.new()
    } else if (ind == max(tmp)){
      plot.new()
    } else {
      real_ind <- real_ind + 1L

      # Offset if different lag
      offset <- (lag - 1) * m

      if (const) {
        offset <- offset + 1
      }

      # i current row, j current column
      i <- rep(1:m, each = m)[real_ind]
      j <- rep(1:m, m)[real_ind]


      curr_main <- ifelse(i == 1, mains[j], "")

      curr_xaxt <- ifelse(i == m, "s", "n") # Only last row gets x axis
      curr_yaxt <- ifelse(j == 1, "s", "n") # Only first column gets y axis

      curr_ylab <- ifelse(j == 1, ylabs[i], "")

      plot(density(x[i, j + offset, ], from = xlim[1], to = xlim[2]),
           xlim = xlim, xaxt = curr_xaxt, yaxt = "n",
           main = "", ylab = curr_ylab, ...)

      par(xpd = NA)
      title(main = curr_main, line = 0.5)
      par(xpd = prev_par$xpd)
    }
  }

  par(prev_par)

}


#' Heatmap of hyperparameters of time-varying coefficient matrix in a TVP-VAR model
#'
#' \code{TV_heatmap} plots a heatmap of  posterior distribution for hyperparameters
#' of a time-varying parameter coefficient matrix in a TVP-VAR model. This is achieved
#' by plotting the median of the posterior of the absolute value for \code{theta_sr} and the
#' median of the posterior for all others. The plot itself is generated by
#' \code{lattice::levelplot}.
#'\code{beta_mean} and \code{theta_sr} will most likely be of interest here.
#'
#'
#' @param x \code{mcmc.var} object
#' @param cuts single integer value, determining the number of cuts for the color palette.
#' The default value is 15.
#' @param cols character string, determining the color palette to be used. The default value is "Purples" for \code{theta_sr}
#' and "RdBu" for all others. See \code{RColorBrewer::brewer.pal.info} for more information.
#' @param max_val numeric value, determining the maximum value for the color palette.
#' If missing, the maximum value is determined by the largest absolute value of the data.
#' @param flipcols logical value, determining whether the color palette should be flipped. The default value is \code{FALSE} for
#' \code{theta_sr} and \code{TRUE} for all others.
#' @param ... further arguments to be passed to \code{lattice::levelplot}.
#'
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' TV_heatmap(res$theta_sr)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
TV_heatmap <- function(x, cuts = 15, cols, max_val, flipcols, ...){

  # Get some attributes
  p <- attr(x, "p")
  const <- attr(x, "const")
  m <- attr(x, "m")
  param_name <- attr(x, "param_name")
  colnames <- attr(x, "colnames")

  # Input checks
  if (!inherits(x, "mcmc.var")) {
    stop("x must be of class mcmc.var")
  }

  if (int_input_bad(cuts) | cuts < 2){
    stop("cuts must be a positive integer larger than 1")
  }

  if (!missing(cols)){
    if (!is.character(cols) | !cols %in% rownames(RColorBrewer::brewer.pal.info)){
      stop("cols must be a character string that is a valid color palette name (see RColorBrewer::brewer.pal.info)")
    }
  }

  if (!missing(max_val)){
    if (numeric_input_bad(max_val)){
      stop("max_val must be a positive numeric value")
    }
  }

  if (!missing(flipcols)){
    if (bool_input_bad(flipcols)){
      stop("flipcols must be a logical value")
    }
  }



  if (param_name == "theta_sr"){
    dat <- abs(x)
    heat <- apply(dat, c(1,2), median)

    if (missing(flipcols)){
      flipcols <- FALSE
    }

  } else {
    dat <- x
    heat <- apply(dat, c(1,2), median)

    if (missing(flipcols)){
      flipcols <- TRUE
    }

  }

  rownames(heat) <- colnames
  col_names <- rep(rownames(heat), p)

  if (const){
    col_names <- c("Int", col_names)
  }
  colnames(heat) <- col_names

  if (missing(max_val)){
    max_val <- max(abs(heat), na.rm = TRUE)
  }


  if (attr(x, "param_name") == "theta_sr"){
    brk <- c(lattice::do.breaks(c(0, max_val), cuts-1), Inf)
    if (missing(cols)) {
      cols <- "Purples"
    }
  } else {
    brk <- c(-Inf, lattice::do.breaks(c(-max_val, max_val), cuts-2), Inf)
    if (missing(cols)){
      cols <- "RdBu"
    }
  }




  if (flipcols == TRUE){
    ind <- 9:1
  } else {
    ind <- 1:9
  }

  col_vals <- brewer.pal(9, cols)[ind]
  col_vals <- colorRampPalette(col_vals)

  seps <- 1:( p - 1) * m
  if (const == TRUE){
    seps <- c(1, seps + 1)
  }
  seps <- seps + 0.5

  p <- levelplot(t(heat[nrow(heat):1,]), col.regions = col_vals, at = brk, ylab = "", xlab = "",
                          panel = function(...){
                            lattice::panel.levelplot(...)
                            lattice::panel.abline(v = seps)
                          }, ...)
  print(p)
  invisible(p)
}


#' Graphical summary of posterior distribution of fitted values for TVP-VAR model
#'
#' \code{plot.shrinkTVPVAR_fit} generates plots visualizing the posterior distribution of fitted values
#' of a model generated by a call to \code{shrinkTVPVAR}.
#'
#' @param x a \code{shrinkTVPVAR_fit} object.
#' @param nplot single integer value, determining the number of plots (i.e. number of equations to visualize at once) to be generated.
#' The default value is 3.
#' @param h_borders numeric vector of length 2, determining the horizontal borders of the plot.
#' The first value is the space between the plot and the left border,
#' the second value is the space between the plot and the right border.
#' Both are fractions of the total width of the plot. The default value is \code{c(0.05, 0.05)}.
#' @param w_borders numeric vector of length 2, determining the vertical borders of the plot.
#' The first value is the space between the plot and the top border,
#' the second value is the space between the plot and the bottom border.
#' Both are fractions of the total height of the plot. The default value is \code{c(0.02, 0.02)}.
#' @param ... further arguments to be passed to \code{\link[graphics]{plot}}.
#'
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' fit <- fitted(res)
#' plot(fit)
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkTVPVAR_fit <- function(x, nplot = 3, h_borders = c(0.05, 0.05), w_borders = c(0.02, 0.02), ...) {

  # Input checks
  if (!inherits(x, "shrinkTVPVAR_fit")) {
    stop("x must be of class shrinkTVPVAR_fit")
  }

  if (int_input_bad(nplot)){
    stop("nplot has to be a single, positive integer")
  }

  if (any(sapply(h_borders, numeric_input_bad_zer)) | sum(h_borders) >= 1 | length(h_borders) > 2 |
      (length(h_borders) == 1 & h_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (any(sapply(w_borders, numeric_input_bad_zer)) | sum(w_borders) >= 1 | length(w_borders) > 2 |
      (length(w_borders) == 1 & w_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  # A bit of a hack to be able to use shrinkTVP's native plotting method
  master_list <- list(x)
  mostattributes(master_list) <- attributes(x)
  attr(master_list[[1]], "type") <- "sample"
  names(master_list) <- "fitted values"
  attr(master_list, "class") <- c("shrinkTVP")

  plot(master_list, pars = c("fitted values"), nplot = nplot, h_borders = h_borders, w_borders = w_borders, ...)
}

#' Graphical summary of posterior predictive density for TVP-VAR-SV model
#'
#' \code{plot.shrinkTVPVAR_forc} generates plots visualizing the posterior predictive density generated by \code{forecast_shrinkTVPVAR}.
#'
#' @param x a \code{shrinkTVPVAR_forc} object.
#' @param nplot single integer value, determining the number of plots (i.e. number of equations to visualize at once) to be generated.
#' The default value is 3.
#' @param h_borders numeric vector of length 2, determining the horizontal borders of the plot.
#' The first value is the space between the plot and the left border,
#' the second value is the space between the plot and the right border.
#' Both are fractions of the total width of the plot. The default value is \code{c(0.05, 0.05)}.
#' @param w_borders numeric vector of length 2, determining the vertical borders of the plot.
#' The first value is the space between the plot and the top border,
#' the second value is the space between the plot and the bottom border.
#' Both are fractions of the total height of the plot. The default value is \code{c(0.02, 0.02)}.
#' @param ... further arguments to be passed to \code{\link[graphics]{plot}}.
#'
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#' forc <- forecast_shrinkTVPVAR(res, n.ahead = 2)
#'
#' plot(forc)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkTVPVAR_forc<- function(x, nplot = 3, h_borders = c(0.05, 0.05), w_borders = c(0.02, 0.02), ...) {

  # Input checks
  if (!inherits(x, "shrinkTVPVAR_forc")) {
    stop("x must be of class shrinkTVPVAR_forc")
  }

  if (int_input_bad(nplot)){
    stop("nplot has to be a single, positive integer")
  }

  if (any(sapply(h_borders, numeric_input_bad_zer)) | sum(h_borders) >= 1 | length(h_borders) > 2 |
      (length(h_borders) == 1 & h_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (any(sapply(w_borders, numeric_input_bad_zer)) | sum(w_borders) >= 1 | length(w_borders) > 2 |
      (length(w_borders) == 1 & w_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  # A bit of a hack to be able to use shrinkTVP's native plotting method
  master_list <- list(x)
  mostattributes(master_list) <- attributes(x)
  attr(master_list[[1]], "type") <- "sample"
  names(master_list) <- "forecast"
  attr(master_list, "class") <- c("shrinkTVP")

  plot(master_list, pars = c("forecast"), nplot = nplot, h_borders = h_borders, w_borders = w_borders, ...)

}

#' Plotting method for \code{mcmc.var} objects
#'
#' @param x \code{mcmc.var} object
#' @param ... further arguments to be passed to \code{\link{TV_heatmap}} function.
#'
#' For further details see \code{\link{TV_heatmap}} function.
#'
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#'
#' plot(res$theta_sr)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.mcmc.var <- function(x, ...) {
  TV_heatmap(x, ...)
}

#' Plotting method for \code{mcmc.tvp.var} objects
#'
#' @param x \code{mcmc.tvp.var} object
#' @param ... further arguments to be passed to \code{\link{state_plotter}} function.
#'
#' For further details see \code{\link{state_plotter}} function.
#'
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#'
#' plot(res$beta)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.mcmc.tvp.var <- function(x, ...) {
  state_plotter(x, ...)
}

#' Plotting method for \code{shrinkTVPVAR} objects
#'
#' @param x \code{shrinkTVPVAR} object
#' @param ... further arguments to be passed to \code{\link{state_plotter}} function.
#'
#' For further details see \code{\link{state_plotter}} function.
#'
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVPVAR(p = 2)
#' data <- sim$data
#'
#' res <- shrinkTVPVAR(data, p = 2)
#'
#' plot(res)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkTVPVAR <- function(x, ...) {
  state_plotter(x$beta, ...)
}

