#ifndef SHRINKDTVP_INT_H
#define SHRINKDTVP_INT_H

#include <RcppArmadillo.h>

// The samplekeeper struct holds all current samples and internal MCMC states
struct samplekeeper_dyn {
  arma::mat beta;
  arma::vec beta_mean;
  arma::vec theta_sr;
  arma::vec tau2;
  arma::vec xi2;
  arma::vec tau2_til;
  arma::vec xi2_til;
  arma::vec kappa2;
  arma::vec lambda2;
  arma::vec kappa2_til;
  arma::vec lambda2_til;
  double kappa2_B;
  double lambda2_B;
  double a_xi;
  double a_tau;
  double c_xi;
  double c_tau;
  double d2;
  double e2;

  arma::vec sv_latent;
  double sv_mu;
  double sv_phi;
  double sv_sigma2;
  double h0;

  arma::vec sigma2;

  bool success;
  std::string fail;

  arma::mat chol_C_N_inv;
  arma::vec m_N;

  arma::mat batches;
  arma::vec curr_sds;
  arma::ivec batch_nrs;
  arma::ivec batch_pos;

  arma::mat psi;
  arma::vec lambda_0;
  arma::mat lambda_p;
  arma::mat kappa_p;
  arma::vec rho_p;
  arma::mat rho_batches;
  arma::vec rho_curr_sds;
  arma::ivec rho_batch_nrs;
  arma::ivec rho_batch_pos;

  bool shrink_inter;
  int inter_column;

  friend std::ostream& operator<<(std::ostream& os, const samplekeeper_dyn& sk) {
    os << "beta: " << sk.beta << "\n";
    os << "beta_mean: " << sk.beta_mean.t() << "\n";
    os << "theta_sr: " << sk.theta_sr.t() << "\n";
    os << "tau2: " << sk.tau2.t() << "\n";
    os << "xi2: " << sk.xi2.t() << "\n";
    os << "tau2_til: " << sk.tau2_til.t() << "\n";
    os << "xi2_til: " << sk.xi2_til.t() << "\n";
    os << "kappa2: " << sk.kappa2.t() << "\n";
    os << "lambda2: " << sk.lambda2.t() << "\n";
    os << "kappa2_til: " << sk.kappa2_til.t() << "\n";
    os << "lambda2_til: " << sk.lambda2_til.t() << "\n";
    os << "kappa2_B: " << sk.kappa2_B << "\n";
    os << "lambda2_B: " << sk.lambda2_B << "\n";
    os << "a_xi: " << sk.a_xi << "\n";
    os << "a_tau: " << sk.a_tau << "\n";
    os << "c_xi: " << sk.c_xi << "\n";
    os << "c_tau: " << sk.c_tau << "\n";
    os << "d2: " << sk.d2 << "\n";
    os << "e2: " << sk.e2 << "\n";
    os << "sv_latent: " << sk.sv_latent.t() << "\n";
    os << "sv_mu: " << sk.sv_mu << "\n";
    os << "sv_phi: " << sk.sv_phi << "\n";
    os << "sv_sigma2: " << sk.sv_sigma2 << "\n";
    os << "h0: " << sk.h0 << "\n";
    os << "sigma2: " << sk.sigma2.t() << "\n";
    os << "success: " << sk.success << "\n";
    os << "fail: " << sk.fail << "\n";
    os << "chol_C_N_inv: " << sk.chol_C_N_inv << "\n";
    os << "m_N: " << sk.m_N.t() << "\n";
    //os << "batches: " << sk.batches << "\n";
    os << "curr_sds: " << sk.curr_sds.t() << "\n";
    os << "batch_nrs: " << sk.batch_nrs.t() << "\n";
    os << "batch_pos: " << sk.batch_pos.t() << "\n";
    os << "psi: " << sk.psi << "\n";
    os << "lambda_0: " << sk.lambda_0.t() << "\n";
    os << "lambda_p: " << sk.lambda_p << "\n";
    os << "kappa_p: " << sk.kappa_p << "\n";
    os << "rho_p: " << sk.rho_p.t() << "\n";
    //os << "rho_batches: " << sk.rho_batches << "\n";
    os << "rho_curr_sds: " << sk.rho_curr_sds.t() << "\n";
    os << "rho_batch_nrs: " << sk.rho_batch_nrs.t() << "\n";
    os << "rho_batch_pos: " << sk.rho_batch_pos.t() << "\n";
    os << "shrink_inter: " << sk.shrink_inter << "\n";
    os << "inter_column: " << sk.inter_column << "\n";
    return os;
  }
};


// The hyperkeeper struct holds all hyperparameters and tuning parameters
struct hyperkeeper_dyn {
  double d1;
  double d2;
  double e1;
  double e2;
  bool learn_lambda2_B;
  bool learn_kappa2_B;
  double lambda2_B;
  double kappa2_B;
  bool learn_a_xi;
  bool learn_a_tau;
  double a_xi;
  double a_tau;
  bool learn_c_xi;
  bool learn_c_tau;
  double c_xi;
  double c_tau;
  bool a_eq_c_xi;
  bool a_eq_c_tau;
  double a_tuning_par_xi;
  double a_tuning_par_tau;
  double c_tuning_par_xi;
  double c_tuning_par_tau;
  double beta_a_xi;
  double beta_a_tau;
  double alpha_a_xi;
  double alpha_a_tau;
  double beta_c_xi;
  double beta_c_tau;
  double alpha_c_xi;
  double alpha_c_tau;
  double Bsigma_sv;
  double a0_sv;
  double b0_sv;
  double bmu;
  double Bmu;
  arma::vec adaptive;
  arma::vec target_rates;
  arma::ivec batch_sizes;
  arma::vec max_adapts;

  // ===== Additional hyperparameters for dynamic TVP =====

  bool iid;               // If true, use iid variance updates; otherwise, non-iid
  arma::vec a_psi;        // Hyperparameter vector for psi sampling (one per regressor)
  arma::vec c_psi;        // Hyperparameter vector for psi sampling (one per regressor)
  double a_rho;           // Hyperparameter for rho sampling
  double b_rho;           // Hyperparameter for rho sampling
  arma::vec alpha_rho;    // Vector of shape parameters for the rho updates (one per regressor)
  arma::vec beta_rho;     // Vector of rate parameters for the rho updates (one per regressor)
  double tuning_par_rho;  // Tuning parameter for the MH step for rho
  bool adaptive_rho;      // Flag to enable adaptive MH for rho sampling
  double target_rate_rho; // Target acceptance rate for the rho MH step
  int batch_size_rho;     // Batch size for adapting rho's proposal distribution
  double max_adapt_rho;   // Maximum adaptation factor for rho's proposal distribution
};


// The dynamic TVP update function (one-step update version)
void shrinkDTVP_int(arma::vec y,
                    arma::mat x,
                    std::string mod_type,
                    hyperkeeper_dyn hyperpara,
                    samplekeeper_dyn& samples);

#endif
