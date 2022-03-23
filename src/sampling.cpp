#include <RcppArmadillo.h>
#include "sampling.h"
#include "myproduct.h"

arma::vec r_norm_vec(arma::vec mean, arma::vec sd) {
  int N = mean.n_rows;
  arma::vec epsilon(N);
  epsilon.randn();
  return mean + sd%epsilon;
}

void sample_mu(arma::mat & mu, arma::vec y, int N, arma::uvec xi, arma::uvec xp, arma::uvec varind, int K, int L, double lambda, double tau){
  for(int k=0; k<K; k++){
    for(int l=0; l<L; l++){
      arma::mat mu_l = mu;
      mu_l.shed_col(l);
      arma::vec mu2x = myprod_skip(N, xi, xp, pow(mu.col(l),2), varind[k], varind[k+1]);
      arma::vec den = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k],varind[k+1]), mu2x);
      arma::vec mux = myprod_skip(N, xi, xp, mu.col(l), varind[k], varind[k+1]);
      arma::vec resid = y - sum(myprod(N, xi, xp, mu_l),1);
      arma::vec muhat = (mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k],varind[k+1]), resid%mux))/(den + tau/lambda);
      mu.col(l).rows(varind[k],varind[k+1]-1) = r_norm_vec(muhat, sqrt(1/(lambda*den + tau)));
    }
  }
}
