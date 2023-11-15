#include <RcppArmadillo.h>
#include "myproduct.h"
#include "up_shape.h"
using namespace Rcpp;
//shape parameters
void up_A(arma::mat & alpha,
          arma::vec & R,
          const arma::mat & loglambda,
          const arma::vec & y,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const double & a){
  arma::mat r =  myprod(y.n_rows, xi, xp, exp(loglambda)); //N,L
  R = sum(r, 1);
  alpha = mysum_t(alpha.n_rows, xi, xp, r.each_col()%(y/R)) + a; //D,L
}