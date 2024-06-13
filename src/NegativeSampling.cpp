#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "NegativeSampling.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec rber(const arma::vec & p){
  int K = p.n_rows;
  arma::vec x = arma::zeros<arma::vec>(K);
  for(int k=0; k<K; k++){
    x(k) = R::rbinom(1, p(k));
  }
  return(x);
}

double NegativeSampling(arma::vec & B,
                        const int & N0, 
                        const arma::vec & probx,
                        const arma::vec & vl){
  for(int m=0; m<N0; m++){
    arma::vec ns =  rber(probx);
    for(int j=0; j<B.n_rows; j++){
      arma::vec nns = ns;
      nns.shed_row(j);
      arma::vec nv = vl;
      nv.shed_row(j);
      B(j) += ns(j) * prod(pow(nv, nns)); 
    }
  }
  return sum(B);
}

//
double NegativeSampling2(arma::vec & B,
                         const int & N0, 
                         const arma::vec & probx,
                         const arma::vec & vl){
  for(int j=0; j<B.n_rows; j++){
    arma::vec nns = probx;
    nns.shed_row(j);
    arma::vec nv = vl;
    nv.shed_row(j);
    B(j) +=  N0 * probx(j) *  prod(1+(nv-1)%nns); 
  }
  return sum(B);
}

/*
arma::vec rbinom_vec(const int N, const arma::vec & p){
  int K = p.n_rows;
  arma::vec x = arma::zeros<arma::vec>(K);
  for(int k=0; k<K; k++){
    x(k) = R::rbinom(N, p(k));
  }
  return(x);
}
*/