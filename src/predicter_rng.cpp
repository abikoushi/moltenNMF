#include <RcppArmadillo.h>
#include "myproduct.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat PoissonGamma_rng(int N,
                           int np,
                           arma::uvec xi,
                           arma::uvec xp,
                           arma::mat alpha,
                           arma::mat beta){
  arma::mat out(N, np);
  int D = alpha.n_rows;
  int L = alpha.n_cols;
  for(int k=0;k<np;k++){
    arma::mat V(D,L);
    for(int i=0;i<D;i++){
      for(int j=0;j<L;j++){
        V(i,j) = arma::randg(arma::distr_param(alpha.at(i,j), beta.at(i,j)));
      }
    }
    arma::vec lam = summyprod(N, xi, xp, V);
    for(int i=0;i<N;i++){
      NumericVector S = Rcpp::rpois(1, arma::as_scalar(lam(i)));
      out(i,k) = S[0];
    }
  }
  return out;
}
 

// [[Rcpp::export]]
arma::mat NegbinGamma_rng(int N,
                            int np,
                            arma::uvec xi,
                            arma::uvec xp,
                            arma::mat alpha,
                            arma::mat beta,
                            double tau){
   arma::mat out(N, np);
   int D = alpha.n_rows;
   int L = alpha.n_cols;
   const double invtau = 1/tau;
   for(int k=0;k<np;k++){
     arma::mat V(D,L);
     for(int i=0;i<D;i++){
       for(int j=0;j<L;j++){
         V(i,j) = arma::randg(arma::distr_param(alpha.at(i,j), beta.at(i,j)));
       }
     }
     arma::vec lam = summyprod(N, xi, xp, V);
     for(int i=0;i<N;i++){
       double r = arma::randg(arma::distr_param(tau, invtau));;
       NumericVector S = Rcpp::rpois(1, r*arma::as_scalar(lam(i)));
       out(i,k) = S[0];
     }
   }
   return out;
 } 