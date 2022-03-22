#include <RcppArmadillo.h>
#include "myproduct.h"

arma::mat mysum_t(int n, arma::uvec xi, arma::uvec xp, arma::mat lam) {
  arma::mat out = arma::zeros<arma::mat>(n,lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(i) += lam.row(xi[j]); 
    }
  }
  return out;
}

arma::mat mysum(int n, arma::uvec xi, arma::uvec xp, arma::mat lam) {
  arma::mat out = arma::zeros<arma::mat>(n,lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(xi[j]) += lam.row(i); 
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat myprod(int n, arma::uvec xi, arma::uvec xp, arma::mat lam) {
  arma::mat out = arma::ones<arma::mat>(n,lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(xi[j]) %= lam.row(i); 
    }
  }
  return out;
}

arma::mat myprod_skip(int n, arma::uvec xi, arma::uvec xp, arma::mat lam, int start, int end) {
  arma::mat out = arma::ones<arma::mat>(n,lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    if(i<start|i>end){
      for(int j=xp[i];j<xp[i+1];j++){
        out.row(xi[j]) %= lam.row(i);
      } 
    }
  }
  return out;
}
