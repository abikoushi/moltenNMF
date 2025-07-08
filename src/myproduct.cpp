#include <RcppArmadillo.h>
#include "myproduct.h"


arma::mat myprod(const int & N,
                 const arma::uvec & xi,
                 const arma::uvec & xp,
                 const arma::mat & lam) {
  arma::mat out = arma::ones<arma::mat>(N, lam.n_cols);
  for(arma::uword i = 0; i < xp.n_rows - 1; i++){
    arma::rowvec li = lam.row(i);
    for(arma::uword j = xp(i); j < xp(i+1); j++){
      out.row(xi(j)) %= li;
    }
  }
  return out;
}

arma::vec myprodvec(const int & n,
                    const arma::uvec & xi, const arma::uvec & xp,
                    const arma::vec & lam) {
  arma::vec out(n, arma::fill::ones);
  for(arma::uword i = 0; i < xp.n_rows - 1; i++){
    double li = lam(i);
    for(arma::uword j = xp(i); j < xp(i+1); j++){
      out(xi(j)) *= li;
    }
  }
  return out;
}


arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam) {
  arma::vec out(n, arma::fill::ones);
  for(int i = start; i < end; i++){
    double li = lam(i);
    for(arma::uword j = xp(i); j < xp(i+1); j++){
      out(xi(j)) *= li;
    }
  }
  return out;
}

arma::mat mysum_t(const int & N,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::mat & lam) {
  arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
  for(arma::uword i = 0; i < xp.n_rows - 1; i++){
    for(arma::uword j = xp(i); j < xp(i + 1); j++){
      out.row(i) += lam.row(xi(j));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::vec summyprod(const int & n, const arma::uvec & xi,
                    const arma::uvec & xp, const arma::mat & lam) {
  arma::vec out = arma::zeros<arma::vec>(n);
  for(arma::uword l = 0; l < lam.n_cols; l++){
    arma::vec pv = arma::ones<arma::vec>(n);
    for(arma::uword i = 0; i < (xp.n_rows - 1); i++){
      for(arma::uword j = xp(i); j < xp(i+1); j++){
        pv(xi(j)) *= lam(i,l); 
      }
    }
    out += pv;
  }
  return out;
}
