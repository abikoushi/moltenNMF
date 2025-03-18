#include "RcppArmadillo.h"
#include "rand.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat rand_init(const arma::mat & alpha, const arma::rowvec & beta){
  arma::mat Z = alpha;
  for (int i=0; i<alpha.n_rows; i++) {
    for (int j=0; j<alpha.n_cols; j++) {
      Z(i,j) = arma::randg(arma::distr_param(alpha(i,j), 1/beta(j)));
    }
  }
  return Z;
}
