#ifndef SHRINKTVP_INT_H
#define SHRINKTVP_INT_H

struct samplekeeper{
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
};


struct hyperkeeper{
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
};

void shrinkTVP_int(arma::vec y,
                   arma::mat x,
                   std::string mod_type,
                   hyperkeeper hyperpara,
                   samplekeeper& samples);

#endif
