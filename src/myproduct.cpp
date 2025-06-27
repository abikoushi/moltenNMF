#include <RcppArmadillo.h>
#include "myproduct.h"


arma::mat myprod(const int & N,
                 const arma::uvec & xi,
                 const arma::uvec & xp,
                 const arma::mat & lam) {
  arma::mat out = arma::ones<arma::mat>(N, lam.n_cols);
  for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      out.row(xi(j)) %= lam.row(i);
    }
  }
  return out;
}

arma::vec myprodvec(const int & n,
                    const arma::uvec & xi, const arma::uvec & xp,
                    const arma::vec & lam) {
  arma::vec out = arma::ones<arma::vec>(n);
  for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      out(xi(j)) *= lam(i);
    }
  }
  return out;
}


arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam) {
  arma::vec out = arma::ones<arma::vec>(n);
  for(int i = start; i < end; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      out(xi(j)) *= lam(i);
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

// arma::mat mysum_t_rv(const int & N,
//                   const arma::uvec & xp0,
//                   const arma::vec & fl) {
//   arma::vec out = arma::zeros<arma::mat>(N);
//   arma::uword M = fl.n_rows;
//   for(arma::uword i = 0; i < xp0.n_rows-1; i++){
//     arma::uword count = xp0(i+1)-xp0(i);
//     if (count == 0) continue;
//     arma::uvec indices = arma::randi<arma::uvec>(count, arma::distr_param(0, M - 1));
//     out.row(i) = sum(fl.rows(indices));
//   }
//   return out;
// }
