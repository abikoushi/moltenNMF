#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_shape.h"
#include "up_rate.h"
#include "ELBO.h"
#include "subset_sp.h"
#include "lr.h"
#include "rand.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List doVB_pois(const arma::vec & y,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const int & D,
               const int & L,
               const int & iter,
               const double & a,
               const double & b,
               arma::mat & V,
               const bool & display_progress){
  int N = y.n_rows;
  //arma::mat V = arma::randg<arma::mat>(D,L);
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  arma::vec ll(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    up_B(N, beta, V, logV, alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, V, logV, R, y, a, b);
    pb.increment();
  }
  ll -= sum(lgamma(y+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}

// [[Rcpp::export]]
List doVB_pois_sp(const int & N,
                  const arma::vec & yv,
                  const arma::uvec & yi,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::uvec & varind,
                  const int & D,
                  const int & L,
                  const int & iter,
                  const double & a,
                  const double & b,
                  arma::mat & V,
                  const bool & display_progress){
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  arma::vec ll(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    up_A_sp(alpha, R, logV, yv, yi, xi, xp, a);
    up_B(N, beta, V, logV, alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois_sp(alpha, beta, V, logV, R, yv, yi, a, b);
    pb.increment();
  }
  ll -= sum(lgamma(yv+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}

/////
//SVB
/////
double up_theta_s_sp(arma::mat &alpha, arma::mat &beta,
                     arma::vec R, arma::mat V, arma::mat logV,
                     const int &N, const arma::uvec &varind,
                     const arma::vec S_yv, const arma::uvec S_yi,
                     const arma::uvec S_xi, const arma::uvec S_xp,
                     const double a, const double b){
  up_A_sp(alpha, R, logV, S_yv, S_yi, S_xi, S_xp, a);
  up_B(N, beta, V, logV, alpha, S_xi, S_xp, varind, b);
  return 0;
}


// [[Rcpp::export]]
List doSVB_pois_sp(const int & N,
                const arma::vec & yv,
                const arma::uvec & yi,
                const arma::uvec & xi,
                const arma::uvec & xp,
                const arma::uvec & varind,
                const int & D,
                const int & L,
                const int & iter,
                const double & a,
                const double & b,
                arma::mat & V,
                const int & bsize,
                const arma::vec & lr_param,
                const std::string & lr_type,
                const bool & display_progress){
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  arma::vec ll = arma::zeros<arma::vec>(iter);
  int N1 =  yi.n_rows;
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  const double NS = ((double) N) / ((double) bsize);
  Progress pb(iter, display_progress);
  for(int epoc = 0; epoc < iter; epoc++){
    arma::umat bags = randpick_c(N, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec S_yi;
      arma::vec S_yv;
      arma::uvec S_xp;
      arma::uvec S_xi;
      filter_y(S_yv, S_yi, yv, yi, bags.col(step));
      subset_spx(S_xi, S_xp, xi, xp, bags.col(step));
      arma::mat alpha_s = alpha;
      arma::mat beta_s = beta;
      arma::vec SR = R.rows(bags.col(step));
      up_theta_s_sp(alpha_s, beta_s,
                    SR, V, logV,
                    bsize, varind,
                    S_yv, S_yi, S_xi, S_xp,
                    a, b);
      ll.row(epoc) += lowerbound_logML_pois_sp(alpha_s, beta_s, V, logV, SR, S_yv, S_yi, a, b);      
      R.rows(bags.col(step)) = SR; 
      beta = rho2 * beta + rho * beta_s * NS;
      alpha = rho2 * alpha + rho * alpha_s * NS;
    }
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}
