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

// [[Rcpp::export]]
List doVB_pois_sp2(const int & N,
                  const arma::vec & yv,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::uvec & varind,
                  const arma::vec & probX0,
                  const double & N0,
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
    up_A(alpha, R, logV, yv, xi, xp, a);
    up_B_sp(N, beta, V, logV, alpha, xi, xp, varind, probX0, N0, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, V, logV, R, yv, a, b);
    pb.increment();
  }
  ll -= sum(lgamma(yv+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}

/*
List doVB_pois_sp2(const int & N,
                   const int & WN,
                  const arma::vec & yv,
                  const arma::uvec & yi,
                  const arma::uvec & xi,
                  const arma::uvec & xp,
                  const arma::vec & weight,
                  const arma::uvec & wxi,
                  const arma::uvec & wxp,
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
    up_B_ww(WN, beta, V, logV, alpha, wxi, wxp, varind, weight, b);
    ll.row(i) = lowerbound_logML_pois_sp(alpha, beta, V, logV, R, yv, yi, a, b);
    pb.increment();
  }
  ll -= sum(lgamma(yv+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}
*/