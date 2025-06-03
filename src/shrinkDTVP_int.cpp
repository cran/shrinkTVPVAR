#include <RcppArmadillo.h>
#include <stochvol.h>
#include <math.h>
#include <shrinkTVP.h>
#include "shrinkDTVP_int.h"       // declares shrinkDTVP_int, hyperkeeper, samplekeeper, etc.
#include "cpp_utilities.h"        // defines res_protector, etc.
using namespace Rcpp;

// This function performs one-step updates for a dynamic TVP model.
void shrinkDTVP_int(arma::vec y,
                    arma::mat x,
                    std::string mod_type,
                    hyperkeeper_dyn hyperpara,
                    samplekeeper_dyn& samples) {

  // Dimensions
  int N = y.n_elem;
  int d = x.n_cols;

  // Override fixed hyperparameter values when learning is turned off
  if (!hyperpara.learn_kappa2_B) samples.kappa2_B = hyperpara.kappa2_B;
  if (!hyperpara.learn_lambda2_B) samples.lambda2_B = hyperpara.lambda2_B;
  if (!hyperpara.learn_a_xi)      samples.a_xi      = hyperpara.a_xi;
  if (!hyperpara.learn_a_tau)     samples.a_tau     = hyperpara.a_tau;
  if (!hyperpara.learn_c_xi)      samples.c_xi      = hyperpara.c_xi;
  if (!hyperpara.learn_c_tau)     samples.c_tau     = hyperpara.c_tau;

  // Update tau2 and xi2 parameters according to model type
  if(mod_type == "triple") {
    shrinkTVP::calc_xi2_tau2(samples.xi2,
                             samples.xi2_til,
                             samples.kappa2_til,
                             samples.kappa2_B,
                             samples.c_xi,
                             samples.a_xi);
    shrinkTVP::calc_xi2_tau2(samples.tau2,
                             samples.tau2_til,
                             samples.lambda2_til,
                             samples.lambda2_B,
                             samples.c_tau,
                             samples.a_tau);
  } else if(mod_type == "ridge") {
    samples.tau2.fill(2.0 / samples.lambda2_B);
    samples.xi2.fill(2.0 / samples.kappa2_B);
  }

  // Define stochastic volatility prior specification
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {
    PriorSpec::Latent0(),
    PriorSpec::Mu(PriorSpec::Normal(hyperpara.bmu, std::sqrt(hyperpara.Bmu))),
    PriorSpec::Phi(PriorSpec::Beta(hyperpara.a0_sv, hyperpara.b0_sv)),
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / hyperpara.Bsigma_sv))
  };
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert = {
    true,
    stochvol::Parameterization::CENTERED,
    1e-8,
    1e-12,
    2,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,
    -1,
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL
  };
  arma::uvec r = arma::uvec(N); r.fill(5);

  arma::mat beta_nc(d, N, arma::fill::zeros);


  // --- Sample time-varying beta (non-centered) ---
  //try {
  // Use the dynamic version of the beta sampler.
  shrinkTVP::sample_beta_McCausland_dyn(beta_nc,
                                        y,
                                        x,
                                        samples.theta_sr,
                                        samples.sigma2,
                                        samples.beta_mean,
                                        samples.psi,
                                        samples.m_N,
                                        samples.chol_C_N_inv);
  //} catch(...) {
  //  beta_nc.fill(std::numeric_limits<double>::quiet_NaN());
  //   if(samples.success) {
  //     samples.fail = "sample_beta_McCausland_dyn failed";
  //     samples.success = false;
  //   }
  // }

  // --- Sample alpha (i.e. beta_mean and theta_sr) ---
  // try {
  shrinkTVP::sample_alpha(samples.beta_mean,
                          samples.theta_sr,
                          y,
                          x,
                          beta_nc,
                          samples.sigma2,
                          samples.tau2,
                          samples.xi2);

  // } catch(...) {
  //   samples.beta_mean.fill(std::numeric_limits<double>::quiet_NaN());
  //   samples.theta_sr.fill(std::numeric_limits<double>::quiet_NaN());
  //   if(samples.success) {
  //     samples.fail = "sample_alpha failed";
  //     samples.success = false;
  //   }
  // }

  // --- Weave into centered parametrization, resample alpha, then revert ---
  shrinkTVP::to_CP(samples.beta, beta_nc, samples.theta_sr, samples.beta_mean);

  //try {
  shrinkTVP::resample_alpha_dyn(samples.beta_mean,
                                samples.theta_sr,
                                samples.beta,
                                beta_nc,
                                samples.psi,
                                samples.xi2,
                                samples.tau2);

  //} catch(...) {
  //  samples.beta_mean.fill(std::numeric_limits<double>::quiet_NaN());
  //  samples.theta_sr.fill(std::numeric_limits<double>::quiet_NaN());
  //  if(samples.success) {
  //    samples.fail = "resample_alpha_dyn failed";
  //    samples.success = false;
  //  }
  //}
  shrinkTVP::to_NCP(beta_nc, samples.beta, samples.theta_sr, samples.beta_mean);




  // --- Sample prior variances (double or triple gamma) ---
  if(mod_type == "double") {
    int placeholder;
    shrinkTVP::sample_DG_TVP(samples.beta_mean,
                             samples.theta_sr,
                             samples.tau2,
                             samples.xi2,
                             samples.lambda2_B,
                             samples.kappa2_B,
                             samples.a_xi,
                             hyperpara.beta_a_xi,
                             hyperpara.alpha_a_xi,
                             samples.a_tau,
                             hyperpara.beta_a_tau,
                             hyperpara.alpha_a_tau,
                             hyperpara.d1,
                             hyperpara.d2,
                             hyperpara.e1,
                             hyperpara.e2,
                             hyperpara.learn_kappa2_B,
                             hyperpara.learn_lambda2_B,
                             hyperpara.learn_a_xi,
                             hyperpara.learn_a_tau,
                             hyperpara.a_tuning_par_xi,
                             hyperpara.a_tuning_par_tau,
                             hyperpara.adaptive,
                             samples.batches,
                             samples.curr_sds,
                             hyperpara.target_rates,
                             hyperpara.max_adapts,
                             samples.batch_nrs,
                             hyperpara.batch_sizes,
                             samples.batch_pos,
                             1, // one-step update mode indicator
                             samples.success,
                             samples.fail,
                             placeholder);
  } else if(mod_type == "triple") {
    int placeholder;
    shrinkTVP::sample_TG_TVP(samples.beta_mean,
                             samples.theta_sr,
                             samples.tau2,
                             samples.xi2,
                             samples.tau2_til,
                             samples.xi2_til,
                             samples.lambda2_til,
                             samples.kappa2_til,
                             samples.lambda2_B,
                             samples.kappa2_B,
                             samples.a_xi,
                             hyperpara.beta_a_xi,
                             hyperpara.alpha_a_xi,
                             samples.a_tau,
                             hyperpara.beta_a_tau,
                             hyperpara.alpha_a_tau,
                             hyperpara.d2,
                             hyperpara.e2,
                             samples.c_xi,
                             samples.c_tau,
                             hyperpara.beta_c_xi,
                             hyperpara.alpha_c_xi,
                             hyperpara.beta_c_tau,
                             hyperpara.alpha_c_tau,
                             hyperpara.learn_kappa2_B,
                             hyperpara.learn_lambda2_B,
                             hyperpara.learn_a_xi,
                             hyperpara.learn_a_tau,
                             hyperpara.learn_c_xi,
                             hyperpara.learn_c_tau,
                             hyperpara.a_tuning_par_xi,
                             hyperpara.a_tuning_par_tau,
                             hyperpara.c_tuning_par_xi,
                             hyperpara.c_tuning_par_tau,
                             hyperpara.a_eq_c_xi,
                             hyperpara.a_eq_c_tau,
                             hyperpara.adaptive,
                             samples.batches,
                             samples.curr_sds,
                             hyperpara.target_rates,
                             hyperpara.max_adapts,
                             samples.batch_nrs,
                             hyperpara.batch_sizes,
                             samples.batch_pos,
                             1, // one-step update mode indicator
                             samples.success,
                             samples.fail,
                             placeholder);
  }

  // --- Sample sigma2 via stochastic volatility ---
  // try {
  // Compute standardized residuals using dynamic beta.
  // (Assumes that columns 1...N of beta_nc correspond to time-specific deviations.)
  arma::vec fitted = x * samples.beta_mean + (x % beta_nc.t()) * samples.theta_sr;
  arma::vec datastand = 2 * arma::log(arma::abs(y - fitted));
  std::for_each(datastand.begin(), datastand.end(), res_protector);

  double sigma = std::sqrt(samples.sv_sigma2);
  arma::vec cur_h = arma::log(samples.sigma2);


  stochvol::update_fast_sv(datastand,
                           samples.sv_mu,
                           samples.sv_phi,
                           sigma,
                           samples.h0,
                           cur_h,
                           r,
                           prior_spec,
                           expert);
  samples.sigma2   = arma::exp(cur_h);
  samples.sv_sigma2 = sigma * sigma;
  std::for_each(samples.sigma2.begin(), samples.sigma2.end(), res_protector);
  // } catch(...) {
  //   samples.sigma2.fill(std::numeric_limits<double>::quiet_NaN());
  //   if(samples.success) {
  //     samples.fail = "sample sigma2 failed";
  //     samples.success = false;
  //   }
  // }

  // --- Variance components: either iid or non-iid updates ---
  if(hyperpara.iid) {
    for (int k = 0; k < d; k++) {
      // Optionally skip the intercept column
      if (!samples.shrink_inter && (k == (samples.inter_column - 1))) continue;
      // try {
      samples.lambda_p(k, 0) = shrinkTVP::sample_lambda_iid(hyperpara.a_psi(k), hyperpara.c_psi(k), samples.psi(k, 0));
      samples.psi(k, 0)      = shrinkTVP::sample_psi_iid(hyperpara.c_psi(k), samples.lambda_p(k, 0), beta_nc(k, 0));
      for (int t = 1; t < N; t++) {
        samples.lambda_p(k, t) = shrinkTVP::sample_lambda_iid(hyperpara.a_psi(k), hyperpara.c_psi(k), samples.psi(k, t));
        samples.psi(k, t)      = shrinkTVP::sample_psi_iid(hyperpara.c_psi(k), samples.lambda_p(k, t),
                    beta_nc(k, t) - beta_nc(k, t-1));
      }
      // } catch(...) {
      //   samples.lambda_p.row(k).fill(std::numeric_limits<double>::quiet_NaN());
      //   samples.psi.row(k).fill(std::numeric_limits<double>::quiet_NaN());
      //   if(samples.success) {
      //     samples.fail = "iid lambda/psi sampling failed";
      //     samples.success = false;
      //   }
      // }
    }
  } else {
    for (int k = 0; k < d; k++) {
      if (!samples.shrink_inter && (k == (samples.inter_column - 1))) continue;
      arma::vec curr_batch;

      if (hyperpara.adaptive_rho) {
        curr_batch = samples.rho_batches.col(k);
      }

      double psi0 = 1.0 / R::rgamma(hyperpara.c_psi(k), 1.0 / samples.lambda_0(k));
      samples.rho_p(k) = shrinkTVP::rho_p_MH_step_marg_oeverything(samples.rho_p(k),
                    samples.psi.row(k).t(),
                    psi0,
                    hyperpara.a_psi(k),
                    hyperpara.c_psi(k),
                    hyperpara.a_rho,
                    hyperpara.b_rho,
                    hyperpara.alpha_rho(k),
                    hyperpara.beta_rho(k),
                    hyperpara.tuning_par_rho,
                    hyperpara.adaptive_rho,
                    curr_batch,
                    samples.rho_curr_sds(k),
                    hyperpara.target_rate_rho,
                    hyperpara.max_adapt_rho,
                    samples.rho_batch_nrs(k),
                    hyperpara.batch_size_rho,
                    samples.rho_batch_pos(k));

      if (hyperpara.adaptive_rho) {
        samples.rho_batches.col(k) = curr_batch;
      }

      // try {
      arma::vec kappa_tmp = samples.kappa_p.row(k).t();
      shrinkTVP::sample_kappa_fast_marg_alternating(kappa_tmp,
                                                    samples.psi.row(k).t(),
                                                    hyperpara.a_psi(k),
                                                    hyperpara.c_psi(k),
                                                    samples.rho_p(k));
      samples.kappa_p.row(k) = kappa_tmp.t();
      // } catch(...) {
      //   samples.kappa_p.row(k).fill(std::numeric_limits<double>::quiet_NaN());
      //   if(samples.success) {
      //     samples.fail = "non-iid kappa_p sampling failed";
      //     samples.success = false;
      //   }
      // }
      // try {
      samples.lambda_0(k) = shrinkTVP::sample_lambda_0(samples.kappa_p(k, 0),
                                                  hyperpara.a_psi(k),
                                                  hyperpara.c_psi(k),
                                                  samples.rho_p(k));
      samples.lambda_p.row(k) = shrinkTVP::sample_lambda(samples.kappa_p.row(k).t(),
                           samples.psi.row(k).t(),
                           hyperpara.a_psi(k),
                           hyperpara.c_psi(k),
                           samples.rho_p(k));
      // } catch(...) {
      //   samples.lambda_p.row(k).fill(std::numeric_limits<double>::quiet_NaN());
      //   if(samples.success) {
      //     samples.fail = "non-iid lambda_p sampling failed";
      //     samples.success = false;
      //   }
      // }
      // try {
      samples.psi.row(k) = shrinkTVP::sample_psi(samples.lambda_p.row(k).t(),
                      beta_nc.row(k).t(),
                      hyperpara.c_psi(k));
      // } catch(...) {
      //   samples.psi.row(k).fill(std::numeric_limits<double>::quiet_NaN());
      //   if(samples.success) {
      //     samples.fail = "non-iid psi sampling failed";
      //     samples.success = false;
      //   }
      // }
    }
  }

  // Weave into centered parameterization again
  // shrinkTVP::to_CP(samples.beta, beta_nc, samples.theta_sr, samples.beta_mean);


  // --- Random sign switch on theta_sr ---
  // for (int i = 0; i < d; i++) {
  //   if (R::runif(0, 1) > 0.5)
  //     samples.theta_sr(i) = -samples.theta_sr(i);
  // }
}
