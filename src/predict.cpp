#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
using namespace Rcpp;

void mapV1(const int & N,
           arma::mat & beta,
           arma::mat & V,
           arma::mat & logV,
           const arma::mat & alpha,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & b,
           int k){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    arma::vec B = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k], varind[k+1]), vl) + b;
    beta.col(l).rows(varind[k], varind[k+1]-1) = B;
    V.col(l).rows(varind[k], varind[k+1]-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
    vl %= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    logV = mat_digamma(alpha) - log(beta);
  }
}