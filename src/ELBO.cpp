#include <RcppArmadillo.h>
#include "ELBO.h"
using namespace Rcpp;

arma::vec elementwise_prod(const arma::vec & yv,
                          const arma::uvec & yi,
                          const arma::vec & R){
  arma::vec out = arma::zeros<arma::vec>(R.n_rows);
  for(int i=0; i<yv.n_rows; i++){
    out(yi(i)) = yv(i)*R(yi(i));
  }
  return out;
}

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

double lowerbound_logML_pois_sp(const arma::mat & alpha,
                             const arma::mat & beta,
                             const arma::mat & lambda,
                             const arma::mat & loglambda,
                             const arma::vec & R,
                             const arma::vec & yv,
                             const arma::uvec & yi,
                             const double & a,
                             const double & b){
  return sum(-R + elementwise_prod(yv, yi, log(R))) - sum(lgamma(yv+1)) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha));
}
