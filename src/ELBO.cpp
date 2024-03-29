#include <RcppArmadillo.h>
#include "ELBO.h"
using namespace Rcpp;

double lowerbound_logML_pois(const arma::mat & alpha,
                             const arma::mat & beta,
                             const arma::mat & lambda,
                             const arma::mat & loglambda,
                             const arma::vec & R,
                             const arma::vec & y,
                             const double & a,
                             const double & b){
  return sum(-R + y%log(R) - lgamma(y+1)) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha));
}
