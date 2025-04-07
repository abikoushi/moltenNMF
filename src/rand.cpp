#include "RcppArmadillo.h"
#include "rand.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

arma::mat rand_init(const arma::mat & alpha, const arma::rowvec & beta){
  arma::mat Z = alpha;
  for (int i=0; i < (int) alpha.n_rows; i++) {
    for (int j=0; j < (int) alpha.n_cols; j++) {
      Z(i,j) = arma::randg(arma::distr_param(alpha(i,j), 1/beta(j)));
    }
  }
  return Z;
}

arma::umat randpick_c(int N1, int b_size){
  if(b_size>N1){
    b_size = N1;    
  }
  arma::uvec rind = arma::randperm(N1);
  int rem =  N1%b_size;
  int sb = N1/b_size;
  arma::umat subind;
  if(rem==0){
    subind = arma::zeros<arma::umat>(b_size, sb);
  }else{
    subind = arma::zeros<arma::umat>(b_size, sb+1);
  }
  if(rem==0){
    for(int i=0; i<sb; i++){
      subind.col(i) = rind.rows(b_size*i, b_size+b_size*i-1);
    }
  }else{
    for(int i=0; i<sb; i++){
      subind.col(i) = rind.rows(b_size*i, b_size+b_size*i-1);
    }
    arma::uvec r_rem  =  rind.rows(0, b_size-rem-1);
    subind.col(sb) = join_cols(rind.rows(b_size*sb, rind.n_rows-1), r_rem);
  }
  return subind;
}