#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
#include "NegativeSampling.h"
// [[Rcpp::depends(RcppArmadillo)]]


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
