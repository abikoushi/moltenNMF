#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate_sp.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

double NegativeSampling(arma::vec & B,
                        double & gamma_l,
                        const int & N0,
                        const arma::vec & probx,
                        const arma::vec & vl){
  double lp = 0;
  for(int j=0; j<B.n_rows; j++){
    double b0 =  N0 * probx(j) *  gamma_l /(1+(vl(j)-1)*probx(j));
    B(j) += b0;
    lp += b0;
  }
  return lp;
}

//rate parameters
double up_B_sp(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          arma::vec & gamma_all,
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
    double gamma_l = gamma_all(l);
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
      arma::vec B = mysum_t(varind(k+1)-varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      //zero-part
      lp += NegativeSampling(B,  gamma_l, N0, probx.rows(varind(k),varind(k+1)-1), V.col(l).rows(varind[k],varind[k+1]-1));
      B += b;
      beta.col(l).rows(varind(k), varind(k+1)-1) = B;
      V.col(l).rows(varind(k), varind(k+1)-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
      //up gamma;
      gamma_l = prod(1+(V.col(l)-1)%probx);
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
    gamma_all(l) = gamma_l;
    logV = mat_digamma(alpha) - log(beta);
  }
  return -lp;
}
