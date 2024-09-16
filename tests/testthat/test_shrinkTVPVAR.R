test_bed <- function(args, m, p) {

  set.seed(123)
  full_dat <- simTVPVAR(N = 20, m = m, p = p, display_progress = FALSE)
  args$y <- full_dat$data

  res <- do.call(shrinkTVPVAR, args)

  expect_s3_class(res, "shrinkTVPVAR")

  # Test prediction methods
  expect_s3_class(fitted(res), "shrinkTVPVAR_fit")
  expect_s3_class(forecast_shrinkTVPVAR(res, 1), "shrinkTVPVAR_forc")

  # Test print method
  expect_visible(res)

  # Test plot methods
  expect_invisible(plot(res))
  param_names <- names(res)[!names(res) %in% c("data", "pred_objs")]
  for (i in param_names) {
    dev.off()
    if (i == "beta_consts") {
      for (j in 1:m) {
        expect_invisible(plot(res[[i]][[j]]))
      }
    } else {
      expect_invisible(plot(res[[i]]))
    }
  }

  forc <- forecast_shrinkTVPVAR(res, 3)
  expect_invisible(plot(forc))

  fit <- fitted(res)
  expect_invisible(plot(fit))
}

mod_type = c("triple", "double", "ridge")
p <- 1:4
m <- 2:4
scenarios <- expand.grid(mod_type, p, m)
names(scenarios) <- c("mod_type", "p", "m")

params <- c(
  "learn_a_xi",
  "learn_a_tau",
  "learn_c_xi",
  "learn_c_tau",
  "a_eq_c_xi",
  "a_eq_c_tau",
  "learn_kappa2_B",
  "learn_lambda2_B",
  "a_xi_adaptive",
  "c_xi_adaptive",
  "a_tau_adaptive",
  "c_tau_adaptive"
)

for(i in length(scenarios)) {

  for (j in params) {

    args <- formals(shrinkTVPVAR)
    args <- args[sapply(args, function(x) x != "")]

    if (!grepl("adaptive", j)) {
      args$TVP_params_beta <- list(tmp = FALSE)
      names(args$TVP_params_beta) <- j

      args$TVP_params_sigma <- list(tmp = FALSE)
      names(args$TVP_params_sigma) <- j
    } else {
      if (j == "a_xi_adaptive") {
        args$TVP_params_beta <- list(adaptive = c(TRUE, FALSE, FALSE, FALSE))
        args$TVP_params_sigma <- list(adaptive = c(TRUE, FALSE, FALSE, FALSE))
      } else if (j == "c_xi_adaptive") {
        args$TVP_params_beta <- list(adaptive = c(FALSE, TRUE, FALSE, FALSE))
        args$TVP_params_sigma <- list(adaptive = c(FALSE, TRUE, FALSE, FALSE))
      } else if (j == "a_tau_adaptive") {
        args$TVP_params_beta <- list(adaptive = c(FALSE, FALSE, TRUE, FALSE))
        args$TVP_params_sigma <- list(adaptive = c(FALSE, FALSE, TRUE, FALSE))
      } else if (j == "c_tau_adaptive") {
        args$TVP_params_beta <- list(adaptive = c(FALSE, FALSE, FALSE, TRUE))
        args$TVP_params_sigma <- list(adaptive = c(FALSE, FALSE, FALSE, TRUE))
      }
    }


    args$p <- scenarios$p[i]
    args$mod_type <- as.character(scenarios$mod_type[i])
    args$display_progress <- FALSE
    args$niter <- 10
    args$nburn <- 0

    test_that(paste0("scenario: ", i,  ", mod_type: ", scenarios$mod_type[i],
                     ", p: ", scenarios$p[i],
                     ", m: ", scenarios$m[i],
                     ", toggled: ", j), {
                       test_bed(args, scenarios$m[i], scenarios$p[i])
                     })

  }
}
