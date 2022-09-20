#include <RcppArmadillo.h>
#include "myproduct.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
List doVB_negbin(arma::vec y,
               arma::uvec xi,
               arma::uvec xp,
               arma::uvec varind,
               int D,
               int L,
               int iter=1000,
               double a=0.5,
               double b=0.001){
  int N = y.n_rows;
  int K = varind.n_rows - 1;
  arma::mat lambda = arma::randg<arma::mat>(D,L);
  arma::mat loglambda = log(lambda);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec ll(iter);
  arma::vec z = arma::ones<arma::vec>(N);
  double tau = 1.0;
  for (int i=0; i<iter; i++) {
    arma::mat r =  myprod(N, xi, xp, exp(loglambda));
    arma::vec den = sum(r,1);
    alpha = mysum_t(D, xi, xp, r.each_col()%(y/den)) + a;
    for(int k=0; k<K; k++){
      arma::mat tmpXl = myprod_skip(N,xi,xp,lambda,varind[k],varind[k+1]);
      tmpXl.each_col() %= z;
      arma::mat B = mysum_t(varind[k+1] - varind[k], xi, xp.rows(varind[k],varind[k+1]), tmpXl) + b;
      beta.rows(varind[k],varind[k+1]-1) = B;
      lambda.rows(varind[k],varind[k+1]-1) = alpha.rows(varind[k],varind[k+1]-1)/B;
    }
    loglambda = mat_digamma(alpha) - log(beta);
    arma::vec ybar = sum(myprod(N,xi,xp,lambda),1);
    arma::vec ahat = (y+tau);
    arma::vec bhat = (ybar+tau);
    arma::vec logz = vec_digamma(ahat) - log(bhat);
    z = ahat/bhat;
    tau += sum(logz - z + log(tau) + 1 - R::digamma(tau))/(N*(1/tau - R::trigamma(tau)));
    ll.row(i) = sum(y +y%log(den) - ybar%z - lgamma(y+1))+
      + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
      - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha))+
      + sum(z%bhat - z*a - ahat%log(bhat) + a*log(a)+lgamma(ahat)-lgamma(a));
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("precision")=tau,
                      Named("ELBO")=ll);
}
