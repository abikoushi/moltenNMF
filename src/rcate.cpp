#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

int rcate0(const int & K){
  int x = 0;
  double U = R::runif(0,1);
  if(U > 1/K){
    for(int k=1; k<K; k++){
      if((k/K)<U & U<=(k+1)/K){
        x = k;
      }
    }
  }
  return(x);
}


int rcate(const vec & p){
  int K = p.n_elem;
  //Rprintf("%d",K);
  vec cump = cumsum(p);
  int x = 0;
  double U = R::runif(0,1);
  if(U > cump(0)){
    for(int k=1; k<K; k++){
      if(cump[k-1]<U & U<=cump[k]){
        x = k;
    }
  }
  }
  return(x);
}


arma::vec rber(const arma::vec & p){
  int K = p.n_rows;
  arma::vec x = arma::zeros<arma::vec>(K);
  for(int k=0; k<K; k++){
  double U = R::runif(0,1);
  if(U < p(k)){
      x(k) = 1.0;
    }
  }
  return(x);
}
