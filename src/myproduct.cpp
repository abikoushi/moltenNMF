#include <RcppArmadillo.h>
#include "myproduct.h"

//////
//Main module for updating shape & rate
//////

// [[Rcpp::export]]
arma::mat myprod(const int & N,
                 const arma::uvec & xi,
                 const arma::uvec & xp,
                 const arma::mat & lam) {
  arma::mat out = arma::ones<arma::mat>(N, lam.n_cols);
  for(int i = 0; i < xp.n_rows - 1; i++){
    for(int j = xp(i); j < xp(i+1); j++){
      out.row(xi(j)) %= lam.row(i);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::vec myprodvec(const int & n,
                    const arma::uvec & xi, const arma::uvec & xp,
                    const arma::vec & lam) {
  arma::vec out = arma::ones<arma::vec>(n);
  for(int i = 0; i < xp.n_rows-1; i++){
    for(int j = xp(i); j < xp(i+1); j++){
      out.row(xi(j)) *= lam(i);
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam) {
  arma::vec out = arma::ones<arma::vec>(n);
  for(int i = start; i < end; i++){
    for(int j = xp(i); j < xp(i+1); j++){
      out.row(xi(j)) *= lam(i);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat mysum_t(const int & N,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::mat & lam) {
  arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
  for(int i = 0; i < xp.n_rows - 1; i++){
    for(int j = xp(i); j < xp(i + 1); j++){
      out.row(i) += lam.row(xi(j));
    }
  }
  return out;
}

//////
//export
//////

// [[Rcpp::export]]
arma::vec summyprod(const int & n, const arma::uvec & xi,
                    const arma::uvec & xp, const arma::mat & lam) {
  arma::vec out = arma::zeros<arma::vec>(n);
  //n_cols=L
  for(int l=0; l<lam.n_cols; l++){
    arma::vec pv = arma::ones<arma::vec>(n);
    for(int i=0; i<xp.n_rows-1; i++){
      for(int j=xp[i];j<xp[i+1];j++){
        pv.row(xi[j]) *= lam(i,l); 
      }
    }
    out += pv;
  }
  return out;
}
