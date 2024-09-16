#include <RcppArmadillo.h>
#include <stochvol.h>
#include <math.h>
#include <shrinkTVP.h>
#include "shrinkTVP_int.h"
#include "cpp_utilities.h"
using namespace Rcpp;



void shrinkTVP_int(arma::vec y,
                   arma::mat x,
                   std::string mod_type,
                   hyperkeeper hyperpara,
                   samplekeeper& samples) {


  // Define some necessary dimensions
  int N = y.n_elem;
  int d = x.n_cols;

  // Override inital values with user specified fixed values
  if (hyperpara.learn_kappa2_B == false){
    samples.kappa2_B = hyperpara.kappa2_B;
  }
  if (hyperpara.learn_lambda2_B == false){
    samples.lambda2_B = hyperpara.lambda2_B;
  }

  if (!hyperpara.learn_a_xi){
    samples.a_xi = hyperpara.a_xi;
  }
  if (!hyperpara.learn_a_tau){
    samples.a_tau = hyperpara.a_tau;
  }

  if (!hyperpara.learn_c_xi){
    samples.c_xi = hyperpara.c_xi;
  }

  if (!hyperpara.learn_c_tau){
    samples.c_tau = hyperpara.c_tau;
  }

  // Initial values and objects holding current values of samples
  arma::mat beta_nc(d, N+1, arma::fill::zeros);
  arma::vec sig2 = arma::exp(samples.sv_latent);

  if (mod_type == "triple") {
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
  } else if (mod_type == "ridge"){
    samples.tau2.fill(2.0/samples.lambda2_B);
    samples.xi2.fill(2.0/samples.kappa2_B);
  }


  // Objects required for stochvol to work
  arma::uvec r(N); r.fill(5);
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {  // prior specification object for the update_*_sv functions
    PriorSpec::Latent0(),  // stationary prior distribution on priorlatent0
    PriorSpec::Mu(PriorSpec::Normal(hyperpara.bmu, std::sqrt(hyperpara.Bmu))),  // normal prior on mu
    PriorSpec::Phi(PriorSpec::Beta(hyperpara.a0_sv, hyperpara.b0_sv)),  // stretched beta prior on phi
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / hyperpara.Bsigma_sv))  // normal(0, Bsigma) prior on sigma
  };  // heavy-tailed, leverage, regression turned off
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert {  // very expert settings for the Kastner, Fruehwirth-Schnatter (2014) sampler
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // sample time varying beta.tilde parameters (NC parametrization)
  try {
    shrinkTVP::sample_beta_McCausland(beta_nc,
                           y,
                           x,
                           samples.theta_sr,
                           samples.sigma2,
                           samples.beta_mean,
                           samples.m_N,
                           samples.chol_C_N_inv);
  } catch (...) {
    beta_nc.fill(nanl(""));
    if (samples.success == true) {
      samples.fail = "sample beta_nc";
      samples.success = false;
    }
  }


  // Sample alpha (theta_sr and beta_mean)
  try {
    shrinkTVP::sample_alpha(samples.beta_mean,
                 samples.theta_sr,
                 y,
                 x,
                 beta_nc,
                 samples.sigma2,
                 samples.tau2,
                 samples.xi2);
  } catch(...){
    samples.beta_mean.fill(nanl(""));
    samples.theta_sr.fill(nanl(""));
    if (samples.success == true){
      samples.fail = "sample alpha";
      samples.success = false;
    }
  }

  // Weave into centered parameterization and resample alpha
  shrinkTVP::to_CP(samples.beta,
        beta_nc,
        samples.theta_sr,
        samples.beta_mean);

  try {
    shrinkTVP::resample_alpha(samples.beta_mean,
                   samples.theta_sr,
                   samples.beta,
                   beta_nc,
                   samples.xi2,
                   samples.tau2);
  } catch(...) {
    samples.beta_mean.fill(nanl(""));
    samples.theta_sr.fill(nanl(""));
    if (samples.success == true) {
      samples.fail = "resample alpha";
      samples.success = false;
    }
  }

  // Weave back into non-centered parameterization
  shrinkTVP::to_NCP(beta_nc,
         samples.beta,
         samples.theta_sr,
         samples.beta_mean);

  // Sample prior variance, differentiating between ridge, double gamma and triple gamma
  if (mod_type == "double") {
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
                  1,
                  samples.success,
                  samples.fail,
                  placeholder);
  } else if (mod_type == "triple") {
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
                  1,
                  samples.success,
                  samples.fail,
                  placeholder);
  }

  // sample sigma2 from homoscedastic or SV case
  try {

    arma::vec datastand = 2 * arma::log(arma::abs(y - x * samples.beta_mean - (x % beta_nc.cols(1,N).t()) * samples.theta_sr));
    std::for_each(datastand.begin(), datastand.end(), res_protector);

    // update_sv needs sigma and not sigma^2
    double sigma = std::sqrt(samples.sv_sigma2);

    arma::vec cur_h = arma::log(samples.sigma2);
    stochvol::update_fast_sv(datastand, samples.sv_mu, samples.sv_phi, sigma, samples.h0, cur_h, r, prior_spec, expert);

    // Write back into sample object
    samples.sigma2 = arma::exp(cur_h);

    // change back to sigma^2
    samples.sv_sigma2 = sigma * sigma;

    std::for_each(samples.sigma2.begin(), samples.sigma2.end(), res_protector);
  } catch(...) {
    samples.sigma2.fill(nanl(""));
    if (samples.success == true) {
      samples.fail = "sample sigma2";
      samples.success = false;
    }
  }


  // Random sign switch
  for (int i = 0; i < d; i++) {
    if (R::runif(0,1) > 0.5) {
      samples.theta_sr(i) = -samples.theta_sr(i);
    }
  }


  // if(mod_type == "triple") {
  //   samples.tau2 = samples.tau2_til;
  //   samples.xi2 = samples.xi2_til;
  // }
}
