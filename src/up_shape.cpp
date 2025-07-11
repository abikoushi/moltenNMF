#include <RcppArmadillo.h>
#include "myproduct.h"
#include "up_shape.h"
using namespace Rcpp;
#include "TimerLogger.h"

// update shape parameters
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

////
//sparse y
////

arma::vec elementwise_div(const arma::vec & yv,
                          const arma::uvec & yi,
                          const arma::vec & R){
  arma::vec out = arma::zeros<arma::vec>(R.n_rows);
  out.rows(yi) = yv/R.rows(yi);
  return out;
}

void up_A_sp(arma::mat & alpha,
             arma::vec & R,
             const arma::mat & loglambda,
             const arma::vec & yv,
             const arma::uvec & yi,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const double & a){
  arma::mat r =  myprod(R.n_rows, xi, xp, exp(loglambda)); //(N, L)
  R = sum(r, 1);
  alpha = mysum_t(alpha.n_rows, xi, xp, r.each_col() % elementwise_div(yv, yi, R)); //D,L
  alpha += a;
}

void up_As_sp(arma::mat & alpha,
             arma::vec & R,
             const arma::mat & loglambda,
             const arma::vec & yv,
             const arma::uvec & yi,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const double & a,
             const double & NS){
  arma::mat r =  myprod(R.n_rows, xi, xp, exp(loglambda)); //(N, L)
  R = sum(r, 1);
  alpha = NS*mysum_t(alpha.n_rows, xi, xp, r.each_col() % elementwise_div(yv, yi, R)) + a;
}

void up_As_sp2(arma::mat & alpha,
              arma::vec & R,
              const arma::mat & loglambda,
              const arma::vec & yv,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const double & a,
              const double & NS){
  //TimerLogger total_timer("upA_total");
  arma::mat r =  myprod(R.n_rows, xi, xp, exp(loglambda)); //(N, L)
  R = sum(r, 1);
  alpha = NS*mysum_t(alpha.n_rows, xi, xp, r.each_col() % (yv/R)) + a;
}
