#include <RcppArmadillo.h>
#include "myproduct.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

arma::mat sweep(int N, arma::uvec xi, arma::uvec xp, arma::vec W){
  int K = W.n_rows;
  arma::mat out = arma::zeros<arma::mat>(N,K);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(xi[j]).col(i) = W.row(i);
    }
  }
  return out;
}

arma::vec sweep2(arma::uvec xi, arma::uvec xp, arma::vec W, arma::vec y){
  int K = W.n_rows;
  arma::vec out = arma::zeros<arma::vec>(K);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(i) += W.row(i)*y.row(xi[j]);
    }
  }
  return out;
}


arma::mat mat_digamma(arma::mat a){
  int K = a.n_rows;
  int L = a.n_cols;
  arma::mat out(K,L);
  for(int k=0;k<K;k++){
    for(int l=0;l<L;l++){
      out(k,l) = R::digamma(a(k,l));
    }
  }
  return out;
}

arma::vec vec_digamma(arma::vec a){
  int K = a.n_rows;
  arma::vec out(K);
  for(int k=0;k<K;k++){
    out[k] = R::digamma(a[k]);
  }
  return out;
}

// [[Rcpp::export]]
List doVB_pois(arma::vec y,
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
  for (int i=0; i<iter; i++) {
    arma::mat U =  myprod(N,xi,xp,exp(loglambda));
    arma::vec den = sum(U,1);
    U.each_col() %= (y/den);
    alpha = mysum_t(D,xi,xp,U) + a;
    for(int k=0; k<K; k++){
      arma::mat tmpXl = myprod_skip(N,xi,xp,lambda,varind[k],varind[k+1]);
      arma::mat B = mysum_t(varind[k+1] - varind[k], xi, xp.rows(varind[k],varind[k+1]), tmpXl) + b;
      beta.rows(varind[k],varind[k+1]-1) = B;
      lambda.rows(varind[k],varind[k+1]-1) = alpha.rows(varind[k],varind[k+1]-1)/B;
    }
    loglambda = mat_digamma(alpha) - log(beta);
    arma::vec ybar = sum(myprod(N,xi,xp,lambda),1);
    ll.row(i) = mean(y +y%log(den) - ybar - lgamma(y+1))+
      + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a))/(D*L) +
      - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha))/(D*L);
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}
