#include <RcppArmadillo.h>
#include "myproduct.h"

//////
//Main module for updating shape & rate
//////

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
      out.row(xi(j)) *= lam(i);
    }
  }
  return out;
}


arma::vec myprodvec_sub(const int & n, const arma::uvec & xi, const arma::uvec & xp,
                        const int & start, const int & end, const arma::vec & lam) {
  arma::vec out = arma::ones<arma::vec>(n);
  for(int i = start; i < end; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      out.row(xi(j)) *= lam(i);
    }
  }
  return out;
}

arma::mat mysum_t(const int & N,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::mat & lam) {
  arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
  for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
    for(int j = xp(i); j < (int) xp(i + 1); j++){
      out.row(i) += lam.row(xi(j));
    }
  }
  return out;
}

// arma::mat mysum_t_ww(const int & N,
//                   const arma::uvec & xi,
//                   const arma::uvec & xp,
//                   const arma::mat & lam,
//                   const arma::vec & weight) {
//   arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
//   for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
//     for(int j = xp(i); j < (int) xp(i + 1); j++){
//       out.row(i) += weight(xi(j)) * lam.row(xi(j));
//     }
//   }
//   return out;
// }

//////
//export
//////

// [[Rcpp::export]]
arma::vec summyprod(const int & n, const arma::uvec & xi,
                    const arma::uvec & xp, const arma::mat & lam) {
  arma::vec out = arma::zeros<arma::vec>(n);
  //n_cols=L
  for(int l = 0; l < (int) lam.n_cols; l++){
    arma::vec pv = arma::ones<arma::vec>(n);
    for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
      for(int j = xp(i); j < (int) xp(i+1); j++){
        pv.row(xi(j)) *= lam(i,l); 
      }
    }
    out += pv;
  }
  return out;
}
