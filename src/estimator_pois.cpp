#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_shape.h"
#include "up_rate.h"
#include "up_rate_sp.h"
#include "ELBO.h"
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
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}

// [[Rcpp::export]]
List doVB_pois_spw(const arma::vec & y,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const int & D,
               const int & N0,
               const arma::vec & probx,
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
  arma::vec gamma_all = arma::zeros<arma::vec>(L);
  for (int l=0;l<L;l++){
    gamma_all.row(l) = prod(1+(V.col(l)-1)%probx);
  }
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    double nlp = up_B_sp(N, beta, V, logV,
                         gamma_all,
                         N0, probx,
                         alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, V, logV, R, y, a, b);
    ll.row(i) += nlp;
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}


// [[Rcpp::export]]
List doVB_pois_s(const arma::vec & y,
                   const arma::uvec & xi,
                   const arma::uvec & xp,
                   const arma::uvec & varind,
                   const int & D,
                   const double & p0,
                   const arma::vec & probx,
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
  arma::vec gamma_all = arma::zeros<arma::vec>(L);
  for (int l=0;l<L;l++){
    gamma_all.row(l) = prod(1+(V.col(l)-1)%probx);
  }
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    double nlp = up_B_s(N, beta, V, logV,
                         gamma_all,
                         p0, probx,
                         alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, V, logV, R, y, a, b);
    ll.row(i) += nlp;
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}
