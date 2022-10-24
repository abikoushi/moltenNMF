#include <RcppArmadillo.h>
#include "myproduct.h"

arma::mat mysum_t(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam) {
  arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i]; j<xp[i+1]; j++){
      out.row(i) += lam.row(xi[j]); 
    }
  }
  return out;
}

arma::mat mysum(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam) {
  arma::mat out = arma::zeros<arma::mat>(N, lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i]; j<xp[i+1]; j++){
      out.row(xi[j]) += lam.row(i); 
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat myprod(const int & N, const arma::uvec & xi, const arma::uvec & xp, const arma::mat & lam) {
  arma::mat out = arma::ones<arma::mat>(N, lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(xi[j]) %= lam.row(i); 
    }
  }
  return out;
}


arma::mat mysum2(int n, arma::uvec xi, arma::uvec xp, arma::mat lam) {
  arma::mat out = arma::zeros<arma::mat>(n, lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i]; j<xp[i+1]; j++){
      out.row(xi[j]) += lam.row(i); 
    }
  }
  return out;
}

arma::mat myprod2(int n, arma::uvec mi, arma::uvec xp, arma::mat lam, arma::vec p) {
  arma::mat out = arma::ones<arma::mat>(n, lam.n_cols);
  for(int i=0; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(mi[j]) %= arma::pow(lam.row(i), p[i]);
    }
  }
  return out;
}


arma::mat myprod_skip(const int & N,
                      const arma::uvec & xi,
                      const arma::uvec & xp,
                      const arma::mat & lam,
                      const int & start,
                      const int & end) {
  arma::mat out = arma::ones<arma::mat>(N, lam.n_cols);
  for(int i=0; i<start; i++){
      for(int j=xp[i];j<xp[i+1];j++){
        out.row(xi[j]) %= lam.row(i);
      }
  }
  for(int i=end; i<xp.n_rows-1; i++){
    for(int j=xp[i];j<xp[i+1];j++){
      out.row(xi[j]) %= lam.row(i);
    }
  }
  return out;
}

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
