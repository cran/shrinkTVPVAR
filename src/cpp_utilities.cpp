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
