#include <RcppArmadillo.h>
#include "myproduct.h"
#include "up_rate2.h"
using namespace Rcpp;

// for negbin and offset factor
void up_B2(const int & N,
           arma::mat & beta,
           arma::mat & V,
           arma::mat & logV,
           const arma::vec & z,
           const arma::mat & alpha,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    vl %= z;
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
      //vl.each_col() %= z;
      arma::vec B = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k], varind[k+1]), vl) + b;
      beta.col(l).rows(varind[k], varind[k+1]-1) = B;
      V.col(l).rows(varind[k], varind[k+1]-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}
