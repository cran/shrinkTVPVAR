#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <float.h>
using namespace Rcpp;

void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 250)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 250) * sign;
  }

  if (std::abs(x) > DBL_MAX * std::pow(10, -250)){
    double sign = std::copysign(1, x);
    x = DBL_MAX * std::pow(10, -250) * sign;
  }

  if (std::isnan(x)){
    throw 1;
  }
}

arma::mat mlag(const arma::mat& X, int p) {
  int Traw = X.n_rows;
  int N = X.n_cols;

  // You really want to make sure Traw > p; otherwise, things will break.
  if (Traw <= p) {
    throw std::invalid_argument("The number of rows must be greater than the lag order p.");
  }

  // Create an output matrix with Traw rows and p*N columns, filled with zeros.
  arma::mat Xlag(Traw, p * N, arma::fill::zeros);

  // For each lag (ii from 1 to p), copy the appropriate submatrix.
  // In R: Xlag[(p+1):Traw, (N*(ii-1)+1):(N*ii)] = X[(p+1-ii):(Traw-ii), 1:N]
  // In C++ (0-indexed): rows [p, Traw-1] receive rows [p-ii, Traw-ii-1] from X.
  for (int ii = 1; ii <= p; ++ii) {
    Xlag.submat(p, N * (ii - 1), Traw - 1, N * ii - 1) =
      X.submat(p - ii, 0, Traw - ii - 1, N - 1);
  }

  return Xlag;
}
