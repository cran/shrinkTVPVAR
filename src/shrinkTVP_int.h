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


  friend std::ostream& operator <<(std::ostream& os, samplekeeper const& a)
  {
    return os << "beta:\n" << a.beta << '\n'
              << "beta_mean:\n" << a.beta_mean << '\n'
              << "theta_sr:\n" << a.theta_sr << '\n'
              << "tau2:\n" << a.tau2 << '\n'
              << "xi2:\n" << a.xi2 << '\n'
              << "tau2_til:\n" << a.tau2_til << '\n'
              << "xi2_til:\n" << a.xi2_til << '\n'
              << "kappa2:\n" << a.kappa2 << '\n'
              << "lambda2:\n" << a.lambda2 << '\n'
              << "kappa2_til:\n" << a.kappa2_til << '\n'
              << "lambda2_til:\n" << a.lambda2_til << '\n'
              << "kappa2_B:\n" << a.kappa2_B << '\n'
              << "lambda2_B:\n" << a.lambda2_B << '\n'
              << "a_xi:\n" << a.a_xi << '\n'
              << "a_tau:\n" << a.a_tau << '\n'
              << "c_xi:\n" << a.c_xi << '\n'
              << "c_tau:\n" << a.c_tau << '\n'
              << "d2:\n" << a.d2 << '\n'
              << "e2:\n" << a.e2 << '\n'
              << "sv_mu:\n" << a.sv_mu << '\n'
              << "sv_phi:\n" << a.sv_phi << '\n'
              << "sv_sigma2:\n" << a.sv_sigma2 << '\n'
              << "h0:\n" << a.h0 << '\n'
              << "sigma2:\n" << a.sigma2 << '\n'
              << "success:\n" << a.success << '\n'
              << "fail:\n" << a.fail << '\n'
              << "chol_C_N_inv:\n" << a.chol_C_N_inv << '\n'
              << "m_N:\n" << a.m_N << '\n'
              << "curr_sds:\n" << a.curr_sds << '\n'
              << "batch_nrs:\n" << a.batch_nrs << '\n'
              << "batch_pos:\n" << a.batch_pos << '\n';
  }
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

  friend std::ostream& operator <<(std::ostream& os, hyperkeeper const& a)
  {
    return os << "d1:\n" << a.d1 << '\n'
              << "d2:\n" << a.d2 << '\n'
              << "e1:\n" << a.e1 << '\n'
              << "e2:\n" << a.e2 << '\n'
              << "learn_lambda2_B:\n" << a.learn_lambda2_B << '\n'
              << "learn_kappa2_B:\n" << a.learn_kappa2_B << '\n'
              << "lambda2_B:\n" << a.lambda2_B << '\n'
              << "kappa2_B:\n" << a.kappa2_B << '\n'
              << "learn_a_xi:\n" << a.learn_a_xi << '\n'
              << "learn_a_tau:\n" << a.learn_a_tau << '\n'
              << "a_xi:\n" << a.a_xi << '\n'
              << "a_tau:\n" << a.a_tau << '\n'
              << "learn_c_xi:\n" << a.learn_c_xi << '\n'
              << "learn_c_tau:\n" << a.learn_c_tau << '\n'
              << "c_xi:\n" << a.c_xi << '\n'
              << "c_tau:\n" << a.c_tau << '\n'
              << "a_eq_c_xi:\n" << a.a_eq_c_xi << '\n'
              << "a_eq_c_tau:\n" << a.a_eq_c_tau << '\n'
              << "a_tuning_par_xi:\n" << a.a_tuning_par_xi << '\n'
              << "a_tuning_par_tau:\n" << a.a_tuning_par_tau << '\n'
              << "c_tuning_par_xi:\n" << a.c_tuning_par_xi << '\n'
              << "c_tuning_par_tau:\n" << a.c_tuning_par_tau << '\n'
              << "beta_a_xi:\n" << a.beta_a_xi << '\n'
              << "beta_a_tau:\n" << a.beta_a_tau << '\n'
              << "alpha_a_xi:\n" << a.alpha_a_xi << '\n'
              << "alpha_a_tau:\n" << a.alpha_a_tau << '\n'
              << "beta_c_xi:\n" << a.beta_c_xi << '\n'
              << "beta_c_tau:\n" << a.beta_c_tau << '\n'
              << "alpha_c_xi:\n" << a.alpha_c_xi << '\n'
              << "alpha_c_tau:\n" << a.alpha_c_tau << '\n'
              << "Bsigma_sv:\n" << a.Bsigma_sv << '\n'
              << "a0_sv:\n" << a.a0_sv << '\n'
              << "b0_sv:\n" << a.b0_sv << '\n'
              << "bmu:\n" << a.bmu << '\n'
              << "Bmu:\n" << a.Bmu << '\n'
              << "adaptive:\n" << a.adaptive << '\n'
              << "target_rates:\n" << a.target_rates << '\n'
              << "batch_sizes:\n" << a.batch_sizes << '\n'
              << "max_adapts:\n" << a.max_adapts << '\n';
  }
};

void shrinkTVP_int(arma::vec y,
                   arma::mat x,
                   std::string mod_type,
                   hyperkeeper hyperpara,
                   samplekeeper& samples);

#endif
