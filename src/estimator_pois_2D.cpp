#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "logexpfuns.h"
#include "rand.h"
#include "up_shape_2D.h"
#include "lr.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

double up_B_2D(const arma::field<arma::mat> & alpha,
               arma::mat & beta,
               arma::field<arma::mat> & V,
               arma::field<arma::mat> & logV,
               const double & b,
               const int & L,
               const arma::uvec & dims, 
               const int & k){
  double lp = 0;
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  for(int l=0; l<L; l++){
    double B0 = sum(V(not_k).col(l));
    lp -= B0;
    beta(k,l) = B0 + b;
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}


double up_theta_2D(arma::field<arma::mat> & alpha,
                   arma::mat & beta,
                   arma::field<arma::mat> & V,
                   arma::field<arma::mat> & logV,
                   const arma::vec & y,
                   const arma::umat & X,
                   const int & L,
                   const arma::uvec & dims, 
                   const double & a,
                   const double & b){
  double lp = 0;
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 0);
  lp += up_B_2D(alpha, beta, V, logV, b, L, dims, 0);
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 1);
  lp += up_B_2D(alpha, beta, V, logV, b, L,dims, 1);
  return lp;
}


// [[Rcpp::export]]
List doVB_pois_2D(arma::field<arma::mat> V,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::uvec & dims,
                  const int & L,
                  const int & iter,
                  const double & a,
                  const double & b, 
                  const bool & display_progress){
  arma::umat X(y.n_rows, 2);
  X.col(0) = rowi;
  X.col(1) = coli;
  //arma::field<arma::mat> V(2);
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  V(0) = arma::randg<arma::mat>(dims(0),L);
  V(1) = arma::randg<arma::mat>(dims(1),L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(dims(0), L);
  alpha(1) = arma::ones<arma::mat>(dims(1), L);
  arma::mat beta = arma::ones<arma::mat>(2, L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    double lp0 = up_theta_2D(alpha, beta, V, logV, y, X, L, dims, a, b);
    lp(i) = lp0 + kld2(alpha, beta, a, b);
    pb.increment();
  }
  lp -= sum(lgamma(y+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}