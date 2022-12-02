#include <RcppArmadillo.h>
#include "myproduct.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

double logsumexp(const arma::vec & x){
  double maxx = max(x);
  double out = maxx + std::log(sum(exp(x-maxx)));
  return out;
}

arma::vec sample_intensity(const int & N,
                           const arma::uvec & xi,
                           const arma::uvec & xp,
                           const arma::mat & alpha,
                           const arma::mat & beta){
  int D = alpha.n_rows;
  int L = alpha.n_cols;
  arma::mat V(D,L);
  for(int i=0;i<D;i++){
    for(int j=0;j<L;j++){
      V(i,j) = arma::randg(arma::distr_param(alpha.at(i,j), 1/beta.at(i,j)));
    }
  }
  return sum(myprod(N, xi, xp, V),1);
}

// [[Rcpp::export]]
arma::mat NegBin_lp(const arma::vec & y,
                    const int & np,
                    const arma::uvec & xi,
                    const arma::uvec & xp,
                    const arma::mat & alpha,
                    const arma::mat & beta,
                    const double & tau){
  int N = y.n_rows;
  const double invtau = 1/tau;
  arma::vec out(N,np);
  for(int k=0;k<np;k++){
    arma::vec lp(N);
    arma::vec lam = sample_intensity(N, xi, xp, alpha, beta);
    for(int i=0;i<N;i++){
      double r = arma::randg(arma::distr_param(tau, invtau));;
      lp.row(i) = R::dpois(y[i], r*arma::as_scalar(lam(i)), true);
    }
    out.col(k) = lp;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat Poisson_lp(const arma::vec & y,
                     const int & np,
                     const arma::uvec & xi,
                     const arma::uvec & xp,
                     const arma::mat & alpha,
                     const arma::mat & beta){
  int N = y.n_rows;
  arma::mat out(N, np);
  for(int k=0;k<np;k++){
    arma::vec lp(N);
    arma::vec lam = sample_intensity(N, xi, xp, alpha, beta);
    for(int i=0;i<N;i++){
      lp.row(i) = R::dpois(y[i], arma::as_scalar(lam(i)), true);
    }
    out.col(k) = lp;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat PoissonGamma_rng(int N,
                           int np,
                           arma::uvec xi,
                           arma::uvec xp,
                           arma::mat alpha,
                           arma::mat beta){
  arma::mat out(N, np);
  for(int k=0;k<np;k++){
    arma::vec lam = sample_intensity(N, xi, xp, alpha, beta);
    for(int i=0;i<N;i++){
      out(i,k) = R::rpois(arma::as_scalar(lam(i)));
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
   const double invtau = 1/tau;
   for(int k=0;k<np;k++){
     arma::vec lam = sample_intensity(N, xi, xp, alpha, beta);
     for(int i=0;i<N;i++){
       double r = arma::randg(arma::distr_param(tau, invtau));;
       out(i,k) =  R::rpois(r*arma::as_scalar(lam(i)));
     }
   }
   return out;
 } 