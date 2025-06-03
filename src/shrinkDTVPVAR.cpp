#include <RcppArmadillo.h>
#include <progress.hpp>
#include <math.h>
#include <shrinkTVP.h>
#include "shrinkDTVP_int.h"   // use the dynamic update version
#include "cpp_utilities.h"
using namespace Rcpp;


//[[Rcpp::export]]
List do_shrinkDTVPVAR(arma::mat y_orig,
                      arma::mat x_orig,
                      std::string mod_type,
                      int niter,
                      int nburn,
                      int nthin,
                      bool display_progress,
                      List TVP_param,
                      List starting_vals)
{

  bool do_dat_G = as<bool>(starting_vals["do_dat_G"]);
  bool inter = as<bool>(starting_vals["const"]);
  int lags = as<int>(starting_vals["lags"]);

  // Get dimensions
  int m = y_orig.n_cols;
  int N = y_orig.n_rows;
  int d = x_orig.n_cols;
  int nsave = std::floor((niter - nburn) / nthin);

  // Starting values and TVP parameters
  arma::cube A_samp = as<arma::cube>(starting_vals["A_st"]);
  arma::cube D_samp = as<arma::cube>(starting_vals["D_st"]);
  List TVP_param_beta = TVP_param["beta"];
  List TVP_param_sigma = TVP_param["sigma"];

  // Bookkeeping objects for each equation
  std::vector<samplekeeper_dyn> samples_beta(m);
  std::vector<samplekeeper_dyn> samples_sigma(m);
  std::vector<hyperkeeper_dyn> hyper_vals_beta(m);
  std::vector<hyperkeeper_dyn> hyper_vals_sigma(m);

  List TVP_samp_beta = starting_vals["beta"];
  List TVP_samp_sigma = starting_vals["sigma"];
  List tmp;
  List tmp_hyper;

  // Transfer values from input lists into our bookkeepers (including dynamic fields)
  for (int j = 0; j < m; j++){
    tmp = TVP_samp_beta["eq" + std::to_string(j)];
    tmp_hyper = TVP_param_beta["eq" + std::to_string(j)];
    samples_beta[j] = {
      /* Basic parameters */
      tmp["beta_st"],                     // beta (centered)
      tmp["beta_mean_st"],                // beta_mean
      tmp["theta_sr_st"],                 // theta_sr
      tmp["tau2_st"],                     // tau2
      tmp["xi2_st"],                      // xi2
      tmp["tau2_st"],                     // tau2_til (initialize equal to tau2)
      tmp["xi2_st"],                      // xi2_til (initialize equal to xi2)
      tmp["kappa2_st"],                   // kappa2
      tmp["lambda2_st"],                  // lambda2
      tmp["kappa2_st"],                   // kappa2_til (initialize equal to kappa2)
      tmp["lambda2_st"],                  // lambda2_til (initialize equal to lambda2)
      tmp["kappa2_B_st"],                 // kappa2_B
      tmp["lambda2_B_st"],                // lambda2_B
      tmp["a_xi_st"],                     // a_xi
      tmp["a_tau_st"],                    // a_tau
      tmp["c_xi_st"],                     // c_xi
      tmp["c_tau_st"],                    // c_tau
      tmp["d2_st"],                       // d2
      tmp["e2_st"],                       // e2

      /* Stochastic volatility parameters */
      tmp["sv_latent_st"],                // sv_latent
      tmp["sv_mu_st"],                    // sv_mu
      tmp["sv_phi_st"],                   // sv_phi
      tmp["sv_sigma2_st"],                // sv_sigma2
      tmp["h0_st"],                       // h0

      /* Observation error variance */
      tmp["sigma2_st"],                   // sigma2

      /* Success flag and failure message */
      true,                               // success
      "NA",                               // fail

      /* Objects for beta sampling routines */
      arma::mat(1, 1, arma::fill::none),  // chol_C_N_inv
      arma::vec(1, arma::fill::none),     // m_N

      /* MCMC adaptation objects */
      arma::mat(arma::max(as<arma::vec>(tmp_hyper["batch_sizes"])), 4, arma::fill::zeros),  // batches
      tmp["tuning_pars"],                 // curr_sds
      arma::ivec(4, arma::fill::ones),    // batch_nrs
      arma::ivec(4, arma::fill::zeros),   // batch_pos

      /* Additional dynamic TVP fields */
      tmp["psi_st"],          // psi: local shrinkage scales for beta increments^
      tmp["lambda_0_st"],
      tmp["lambda_p_st"],     // lambda_p: additional shrinkage parameter
      tmp["kappa_p_st"],      // kappa_p: extra shrinkage parameter for non-iid cases
      tmp["rho_p_st"],        // rho_p: additional parameter for non-iid updates
      tmp["rho_batches_st"],  // rho_batches: batches for the rho updates
      tmp["rho_curr_sds_st"], // rho_curr_sds: current standard deviations for rho updates
      tmp["rho_batch_nrs_st"],// rho_batch_nrs: batch numbers for rho updates
      tmp["rho_batch_pos_st"],// rho_batch_pos: batch positions for rho updates

      /* Intercept controls */
      tmp["shrink_inter_st"], // shrink_inter: whether to shrink the intercept
      tmp["inter_column_st"]  // inter_column: 1-indexed intercept column
    };

    tmp = TVP_samp_sigma["eq" + std::to_string(j)];
    tmp_hyper = TVP_param_sigma["eq" + std::to_string(j)];
   samples_sigma[j] = {
     tmp["beta_st"],
     tmp["beta_mean_st"],
     tmp["theta_sr_st"],
     tmp["tau2_st"],
     tmp["xi2_st"],
     tmp["tau2_st"],
     tmp["xi2_st"],
     tmp["kappa2_st"],
     tmp["lambda2_st"],
     tmp["kappa2_st"],
     tmp["lambda2_st"],
     tmp["kappa2_B_st"],
     tmp["lambda2_B_st"],
     tmp["a_xi_st"],
     tmp["a_tau_st"],
     tmp["c_xi_st"],
     tmp["c_tau_st"],
     tmp["d2_st"],
     tmp["e2_st"],
     tmp["sv_latent_st"],
     tmp["sv_mu_st"],
     tmp["sv_phi_st"],
     tmp["sv_sigma2_st"],
     tmp["h0_st"],
     tmp["sigma2_st"],
     true,
     "NA",
     arma::mat(1, 1, arma::fill::none),
     arma::vec(1, arma::fill::none),
     arma::mat(arma::max(as<arma::vec>(tmp_hyper["batch_sizes"])), 4, arma::fill::zeros),
     tmp["tuning_pars"],
     arma::ivec(4, arma::fill::ones),
     arma::ivec(4, arma::fill::zeros),
     tmp["psi_st"],          // psi
     tmp["lambda_0_st"],
     tmp["lambda_p_st"],     // lambda_p
     tmp["kappa_p_st"],      // kappa_p
     tmp["rho_p_st"],        // rho_p
     tmp["rho_batches_st"],  // rho_batches
     tmp["rho_curr_sds_st"], // rho_curr_sds
     tmp["rho_batch_nrs_st"],// rho_batch_nrs
     tmp["rho_batch_pos_st"],// rho_batch_pos
     tmp["shrink_inter_st"], // shrink_inter
     tmp["inter_column_st"]  // inter_column
   };

    // Transfer hyperparameters (they remain as in the non-dynamic VAR)
    tmp = TVP_param_beta["eq" + std::to_string(j)];
   hyper_vals_beta[j] = {
     tmp["d1"],
     tmp["d2"],
     tmp["e1"],
     tmp["e2"],
     tmp["learn_lambda2_B"],
     tmp["learn_kappa2_B"],
     tmp["lambda2_B"],
     tmp["kappa2_B"],
     tmp["learn_a_xi"],
     tmp["learn_a_tau"],
     tmp["a_xi"],
     tmp["a_tau"],
     tmp["learn_c_xi"],
     tmp["learn_c_tau"],
     tmp["c_xi"],
     tmp["c_tau"],
     tmp["a_eq_c_xi"],
     tmp["a_eq_c_tau"],
     tmp["a_tuning_par_xi"],
     tmp["a_tuning_par_tau"],
     tmp["c_tuning_par_xi"],
     tmp["c_tuning_par_tau"],
     tmp["beta_a_xi"],
     tmp["beta_a_tau"],
     tmp["alpha_a_xi"],
     tmp["alpha_a_tau"],
     tmp["beta_c_xi"],
     tmp["beta_c_tau"],
     tmp["alpha_c_xi"],
     tmp["alpha_c_tau"],
     tmp["Bsigma_sv"],
     tmp["a0_sv"],
     tmp["b0_sv"],
     tmp["bmu"],
     tmp["Bmu"],
     tmp["adaptive"],
     tmp["target_rates"],
     tmp["batch_sizes"],
     tmp["max_adapts"],

     // New dynamic components
     tmp["iid"],
     tmp["a_psi"],
     tmp["c_psi"],
     tmp["a_rho"],
     tmp["b_rho"],
     tmp["alpha_rho"],
     tmp["beta_rho"],
     tmp["tuning_par_rho"],
     tmp["adaptive_rho"],
     tmp["target_rate_rho"],
     tmp["batch_size_rho"],
     tmp["max_adapt_rho"]
   };

    tmp = TVP_param_sigma["eq" + std::to_string(j)];
   hyper_vals_sigma[j] = {
     tmp["d1"],
     tmp["d2"],
     tmp["e1"],
     tmp["e2"],
     tmp["learn_lambda2_B"],
     tmp["learn_kappa2_B"],
     tmp["lambda2_B"],
     tmp["kappa2_B"],
     tmp["learn_a_xi"],
     tmp["learn_a_tau"],
     tmp["a_xi"],
     tmp["a_tau"],
     tmp["learn_c_xi"],
     tmp["learn_c_tau"],
     tmp["c_xi"],
     tmp["c_tau"],
     tmp["a_eq_c_xi"],
     tmp["a_eq_c_tau"],
     tmp["a_tuning_par_xi"],
     tmp["a_tuning_par_tau"],
     tmp["c_tuning_par_xi"],
     tmp["c_tuning_par_tau"],
     tmp["beta_a_xi"],
     tmp["beta_a_tau"],
     tmp["alpha_a_xi"],
     tmp["alpha_a_tau"],
     tmp["beta_c_xi"],
     tmp["beta_c_tau"],
     tmp["alpha_c_xi"],
     tmp["alpha_c_tau"],
     tmp["Bsigma_sv"],
     tmp["a0_sv"],
     tmp["b0_sv"],
     tmp["bmu"],
     tmp["Bmu"],
     tmp["adaptive"],
     tmp["target_rates"],
     tmp["batch_sizes"],
     tmp["max_adapts"],

     // New dynamic components
     tmp["iid"],
     tmp["a_psi"],
     tmp["c_psi"],
     tmp["a_rho"],
     tmp["b_rho"],
     tmp["alpha_rho"],
     tmp["beta_rho"],
     tmp["tuning_par_rho"],
     tmp["adaptive_rho"],
     tmp["target_rate_rho"],
     tmp["batch_size_rho"],
     tmp["max_adapt_rho"]
   };
  }


  // Storage objects (identical to the non-dynamic VAR)
  arma::cube beta_store(m * d, N, nsave);
  arma::cube beta_const_store(m, N, nsave);
  arma::cube beta_mean_store(m, d, nsave);
  arma::cube xi2_store(m, d, nsave);
  arma::cube c_xi_store(m, 1, nsave);
  arma::cube kappa2_store(m, d, nsave);
  arma::cube kappa2_B_store(m, 1, nsave);
  arma::cube a_xi_store(m, 1, nsave);
  arma::cube tau2_store(m, d, nsave);
  arma::cube c_tau_store(m, 1, nsave);
  arma::cube lambda2_store(m, d, nsave);
  arma::cube lambda2_B_store(m, 1, nsave);
  arma::cube a_tau_store(m, 1, nsave);
  arma::cube theta_sr_store(m, d, nsave);
  arma::cube SIGMA_store(m * N, m, nsave);

  arma::cube final_A_store(m, m, nsave);
  arma::mat final_D_store(nsave, m);
  arma::cube theta_sr_SIGMA_store(m, m, nsave);
  arma::mat sv_phi_store(nsave, m);
  arma::mat sv_mu_store(nsave, m);
  arma::mat sv_sigma2_store(nsave, m);

  // New storage objects for dynamic case
  arma::cube rho_p_store(m, d, nsave);
  arma::cube psi_store(m * d, N, nsave);
  arma::cube psi_const_store(m, N, nsave);
  arma::cube final_lambda_store(m, d, nsave);

  arma::cube final_lambda_SIGMA_store(m, m, nsave);
  arma::cube rho_p_SIGMA_store(m, m, nsave);

  int stor_pos;

  // Progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // Auxiliary objects
  arma::mat resids(N, m, arma::fill::ones);
  arma::vec aug_check(N, arma::fill::zeros);
  arma::vec y_check;
  arma::vec y_star;
  arma::mat curr_beta;
  arma::mat curr_A;
  List curr_beta_eq_vals;
  List curr_beta_hyp_vals;
  List curr_sigma_eq_vals;
  List curr_sigma_hyp_vals;

  int fail_iter = 0;
  std::string fail;
  int fail_eq = 999;
  std::string fail_comp;
  bool succesful = true;
  std::string int_fail;
  List internals;

  arma::mat y = y_orig;
  arma::mat x = x_orig;

  // Begin Gibbs sampling loop
  for (int irep = 0; irep < niter; irep++)
  {
    // Equation-by-equation updates for beta then sigma
    for (int eq = 0; eq < m; eq++)
    {
      // Build augmented regressor using previous equations
      aug_check.fill(0);
      if (eq > 0)
      {
        for (int j = 0; j < eq; j++)
        {
          for (int t = 0; t < N; t++)
          {
            aug_check(t) += (A_samp.slice(t))(eq, j) * (resids.col(j))(t);
          }
        }
      }

      y_check = y.col(eq) - aug_check;
      // Use the dynamic update function for beta
      shrinkDTVP_int(y_check,
                     x,
                     mod_type,
                     hyper_vals_beta[eq],
                     samples_beta[eq]);
      if (!samples_beta[eq].success){
        succesful = false;
        fail_iter = irep + 1;
        fail = samples_beta[eq].fail;
        fail_eq = eq;
        fail_comp = "beta";
        break;
      }

      // Update sigma parameters for this equation
      curr_beta = (samples_beta[eq].beta);


      if (eq == 0)
      {
        resids.col(0) = y.col(0) - arma::sum(x % curr_beta.t(), 1);
        for (int t = 0; t < N; t++){
          (D_samp.slice(t))(eq, eq) = samples_beta[eq].sigma2(t);
        }
      } else {
        y_star = y.col(eq) - arma::sum(x % curr_beta.t(), 1);
        arma::mat x_star = resids.cols(0, eq - 1);
        shrinkDTVP_int(y_star,
                       x_star,
                       mod_type,
                       hyper_vals_sigma[eq],
                       samples_sigma[eq]);

        if (!samples_sigma[eq].success){
          succesful = false;
          fail_iter = irep + 1;
          fail = samples_sigma[eq].fail;
          fail_eq = eq;
          fail_comp = "sigma";
          break;
        }

        curr_A = (samples_sigma[eq].beta);
        resids.col(eq) = y_star - arma::sum(x_star % curr_A.t(), 1);
        for (int t = 0; t < N; t++){
          for (int j = 0; j < eq; j++){
            (A_samp.slice(t))(eq, j) = curr_A(j, t);
          }
          (D_samp.slice(t))(eq, eq) = samples_sigma[eq].sigma2(t);
        }
        samples_beta[eq].sigma2 = samples_sigma[eq].sigma2;
      }


      // Store outputs every nthin iteration after burn-in
      if ((irep % nthin == 0) && (irep >= nburn))
      {

        stor_pos = irep - nburn;

        (beta_mean_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].beta_mean);
        (theta_sr_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].theta_sr);
        (xi2_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].xi2);
        (c_xi_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].c_xi;
        (kappa2_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].kappa2);
        (kappa2_B_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].kappa2_B;
        (a_xi_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].a_xi;
        (tau2_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].tau2);
        (c_tau_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].c_tau;
        (lambda2_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].lambda2);
        (lambda2_B_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].lambda2_B;
        (a_tau_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].a_tau;
        (beta_store.slice(stor_pos/nthin)).rows(d * eq, d - 1 + (d * eq)) = curr_beta;
        (rho_p_store.slice(stor_pos/nthin)).row(eq) = arma::trans(samples_beta[eq].rho_p);
        (psi_store.slice(stor_pos/nthin)).rows(d * eq, d - 1 + (d * eq)) = samples_beta[eq].psi;

        if (inter) {
          (beta_const_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].beta.submat(0, 0, 0, N-1);
          (psi_const_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].psi.submat(0, 0, 0, N-1);
        }

        if (eq == 0) {
          sv_phi_store(stor_pos/nthin, eq) = samples_beta[eq].sv_phi;
          sv_mu_store(stor_pos/nthin, eq) = samples_beta[eq].sv_mu;
          sv_sigma2_store(stor_pos/nthin, eq) = samples_beta[eq].sv_sigma2;
        } else {
          sv_phi_store(stor_pos/nthin, eq) = samples_sigma[eq].sv_phi;
          sv_mu_store(stor_pos/nthin, eq) = samples_sigma[eq].sv_mu;
          sv_sigma2_store(stor_pos/nthin, eq) = samples_sigma[eq].sv_sigma2;
        }

        (final_lambda_store.slice(stor_pos/nthin)).row(eq) = samples_beta[eq].lambda_p.col(N-1).t();

        for (int j = 0; j < eq; j++) {
          (theta_sr_SIGMA_store.slice(stor_pos/nthin))(eq, j) = samples_sigma[eq].theta_sr[j];
          (final_lambda_SIGMA_store.slice(stor_pos/nthin))(eq, j) = (samples_sigma[eq].lambda_p.col(N-1))[j];
          (rho_p_SIGMA_store.slice(stor_pos/nthin))(eq, j) = samples_sigma[eq].rho_p[j];
        }

        if (eq == (m - 1)) {
          for (int t = 0; t < N; t++) {
            (SIGMA_store.slice(stor_pos/nthin)).rows(m * t, (m - 1) + m * t) =
              A_samp.slice(t) * D_samp.slice(t) * (A_samp.slice(t)).t();
            if (t == (N - 1)) {
              final_A_store.slice(stor_pos/nthin) = A_samp.slice(t);
              final_D_store.row(stor_pos/nthin) = ((D_samp.slice(t)).diag()).t();
            }
          }
        }
      }
    }

    if (do_dat_G) {
      arma::mat y_aug(N + lags, m, arma::fill::zeros);
      int offset = 0;
      for (int t = 0; t < lags; t++){
        arma::vec fill = Rcpp::rnorm(m, 0, 1);
        y_aug.row(t) = fill.t();
      }
      for (int t = lags; t < (N + lags); t++){
        if (inter){
          arma::rowvec curr_inter(m);
          for (int eq = 0; eq < m; eq++){
            curr_inter(eq) = samples_beta[eq].beta(0, t - lags);
          }
          y_aug.row(t) += curr_inter;
          offset = 1;
        }
        for (int lag = 1; lag < (lags + 1); lag++){
          arma::mat curr_PHI(m, m, arma::fill::zeros);
          for (int eq = 0; eq < m; eq++){

            curr_PHI.row(eq) = (samples_beta[eq].beta.submat(((lag - 1) * m) + offset,
                                t - lags,
                                lag * m - 1 + offset,
                                t - lags)).t();

          }
          y_aug.row(t) += (curr_PHI * (y_aug.row(t - lag)).t()).t();
        }
        arma::vec fill = Rcpp::rnorm(m, 0, 1);
        y_aug.row(t) += (A_samp.slice(t - lags) * arma::sqrt(D_samp.slice(t - lags)) * fill).t();
      }
      arma::mat y_lagged = mlag(y_aug, lags);
      if (inter){
        arma::vec intercept(y_lagged.n_rows, arma::fill::ones);
        y_lagged = arma::join_rows(intercept, y_lagged);
      }
      y = y_aug.rows(lags, y_aug.n_rows - 1);
      x = y_lagged.rows(lags, y_lagged.n_rows - 1);
    }

    if (arma::any(prog_rep_points == irep) && display_progress)
      p.increment();

    if (irep % 500 == 0)
      Rcpp::checkUserInterrupt();

    if (!succesful)
      break;
  }

  return List::create(_["beta"] = beta_store,
                      _["beta_const"] = beta_const_store,
                      _["beta_mean"] = beta_mean_store,
                      _["theta_sr"] = theta_sr_store,
                      _["xi2"] = xi2_store,
                      _["c_xi"] = c_xi_store,
                      _["kappa2"]  = kappa2_store,
                      _["kappa2_B"] = kappa2_B_store,
                      _["a_xi"] = a_xi_store,
                      _["tau2"] = tau2_store,
                      _["c_tau"] = c_tau_store,
                      _["lambda2"] = lambda2_store,
                      _["lambda2_B"] = lambda2_B_store,
                      _["a_tau"] = a_tau_store,
                      _["Sigma"] = SIGMA_store,
                      _["psi"] = psi_store,
                      _["psi_const"] = psi_const_store,
                      _["rho_p"] = rho_p_store,
                      _["pred_objs"] = List::create(
                        _["final_A"] = final_A_store,
                        _["final_D"] = final_D_store,
                        _["theta_sr_SIGMA"] = theta_sr_SIGMA_store,
                        _["sv_sigma2"] = sv_sigma2_store,
                        _["sv_phi"] = sv_phi_store,
                        _["sv_mu"] = sv_mu_store),
                        _["final_lambda"] = final_lambda_store,
                        _["final_lambda_SIGMA"] = final_lambda_SIGMA_store,
                        _["rho_p_SIGMA"] = rho_p_SIGMA_store,
                        _["success_vals"] = List::create(
                          _["success"] = succesful,
                          _["fail_iter"] = fail_iter,
                          _["fail"] = fail,
                          _["fail_eq"] = fail_eq,
                          _["fail_comp"] = fail_comp));
}
