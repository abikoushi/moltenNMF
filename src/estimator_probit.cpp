#include <RcppArmadillo.h>
#include <math.h>
#include "myproduct.h"
#include "sampling.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

static const double t4 = 0.45;

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a) {
  double ainv = 1.0 / a;
  double x;
  double rho;
  do {
    x = R::rexp(ainv) + a; /* rexp works with 1/lambda */
    rho = exp(-0.5 * pow((x - a), 2));
  } while (R::runif(0, 1) > rho);
  return x;
}

/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a) {
  double x = -DBL_MAX;
  while (x < a) {
    x = R::rnorm(0, 1);
  }
  return x;
}

double rhalf_norm(double mean, double sd) {
  const double alpha = (- mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha);
  } else {
    return mean + sd * ers_a_inf(alpha);
  }
}

arma::vec rhalf_norm_vec(arma::vec y, arma::vec mean) {
  int N = mean.n_rows;
  arma::vec out(N);
  for(int n=0;n<N;n++){
  	if(y[n]==1){
  		out.row(n) = rhalf_norm(mean[n], 1);
  	}else{
  		out.row(n) = rhalf_norm(-mean[n], 1);  			
  	}
  }
  return out;
}

// [[Rcpp::export]]
List doGibbs_probit(arma::ivec y,
             arma::uvec xi,
             arma::uvec xp,
             arma::uvec varind,
             int D,
             int L,
             int iter=1000,
             double lambda=1,
             double tau=1){
  int N = y.n_rows;
  int K = varind.n_rows - 1;
  arma::mat mu = arma::randn<arma::mat>(D,L);
  arma::mat loglik = arma::zeros<arma::mat>(N,iter);
  arma::cube mu_s(D,L,iter);
  arma::vec f =  sum(myprod(N, xi, xp, mu),1);
  arma::vec Lam(N);
  Lam.fill(lambda);
  for(int i=0; i<iter; i++){
  	arma::vec z =  rhalf_norm_vec(f, Lam);
    sample_mu(mu, z, N, xi, xp, varind, K, L, lambda, tau);
    f =  sum(myprod(N, xi, xp, mu),1);
    loglik.col(i) = y% log(normcdf(f)) + (1-y)%log(normcdf(-f));
    mu_s.slice(i) = mu;
  }
  return List::create(Named("mu") = mu_s, Named("loglik") = loglik);
}

// [[Rcpp::export]]
List doMCVB_probit(arma::ivec y,
                    arma::uvec xi,
                    arma::uvec xp,
                    arma::uvec varind,
                    int D,
                    int L,
                    int iter=1000,
                    double lambda=1,
                    double tau=1){
  int N = y.n_rows;
  int K = varind.n_rows - 1;
  arma::mat m = arma::randn<arma::mat>(D,L);
  arma::mat s2 = arma::ones<arma::mat>(D,L);
  arma::vec loglik = arma::zeros<arma::vec>(iter);
  arma::vec f =  sum(myprod(N, xi, xp, m),1);
  arma::vec Lam(N);
  Lam.fill(lambda);
  for(int i=0; i<iter; i++){
    arma::vec z =  rhalf_norm_vec(f, Lam);
    up_ms(m, s2, z, N, xi, xp, varind, K, L, lambda, tau);
    f =  sum(myprod(N, xi, xp, m),1);
    loglik.row(i) = mean(y% log(normcdf(f)) + (1-y)%log(normcdf(-f)));
  }
  return List::create(Named("mean") = m, Named("sd") = sqrt(s2), Named("loglik") = loglik);
}
