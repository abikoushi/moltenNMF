#include <RcppArmadillo.h>
#include "myproduct.h"
#include "up_shape.h"
using namespace Rcpp;

arma::vec elementwise_div(const arma::vec & yv,
                          const arma::uvec & yi,
                          const arma::vec & R){
  arma::vec out = arma::zeros<arma::vec>(R.n_rows);
  // for(int i=0; i<yv.n_rows; i++){
  //   out(yi(i)) = yv(i)/R(yi(i));
  // }
  out.rows(yi) = yv/R.rows(yi);
  return out;
}

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

//shape parameters
void up_A_sp(arma::mat & alpha,
             arma::vec & R,
             const arma::mat & loglambda,
             const arma::vec & yv,
             const arma::uvec & yi,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const double & a){
  arma::mat r =  myprod(R.n_rows, xi, xp, exp(loglambda)); //N,L
  R = sum(r, 1);
  alpha = mysum_t(alpha.n_rows, xi, xp, r.each_col()%elementwise_div(yv,yi,R)); //D,L
  alpha += a;
}
