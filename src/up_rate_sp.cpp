#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate_sp.h"
#include "NegativeSampling.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//rate parameters
double up_B_sp(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const int & N0,
          const arma::vec & probx,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  double lp = 0.0; 
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
      arma::vec B = mysum_t(varind(k+1)-varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      //sampling
      lp += NegativeSampling2(B, N0, probx, V.col(l));
      B += b;
      beta.col(l).rows(varind[k], varind[k+1]-1) = B;
      V.col(l).rows(varind[k], varind[k+1]-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
  return -lp;
}
