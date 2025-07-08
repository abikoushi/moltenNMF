#include <RcppArmadillo.h>
#include <random>
#include "geometricsampling.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//geometric sampling for SVB
arma::vec geomprod_all(const int & M,
                       const arma::vec & Vl,
                       const arma::uvec & varind,
                       arma::umat & U){
  arma::vec fl(M, arma::fill::ones); // sampling X with size M
  for(arma::uword i = 0; i < (varind.n_rows - 1); i++){
    arma::uvec indices;
    if (varind(i) == varind(i+1) - 1) {
      indices = arma::uvec(M, arma::fill::value(varind(i)));
    } else {
      const int start = varind(i);
      const int end = (varind(i+1) - 1);
      indices = arma::randi<arma::uvec>(M, arma::distr_param(start, end));
    }
    fl %= Vl(indices);
    U.col(i) = indices;
  }
  return fl;
}

arma::vec geomsum_k(const int & D,
                    const int & M, 
                    const double & MR,
                    const arma::vec & fl,
                    const arma::vec & Vl,
                    const arma::uvec & varind,
                    const arma::uword & k,
                    const arma::umat & U){
  arma::vec out(D, arma::fill::zeros);
  if (M == 0) {
    return out;
  }
  arma::uvec out_indices = U.col(k) - varind(k);
  for(int j=0; j < M; j++){
    arma::uword oj = out_indices(j);
    out(oj) += MR * sum(fl/Vl(oj));
  }
  return out;
}
