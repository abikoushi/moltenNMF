#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_shape.h"
#include "up_rate2.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//precision parameter
void up_tau(double & tau,
            arma::vec & z,
            const arma::vec & y,
            const arma::mat & V,
            const arma::uvec & xi,
            const arma::uvec & xp){
  arma::vec ybar = sum(myprod(z.n_rows, xi, xp, V), 1);
  arma::vec ahat = (y+tau);
  arma::vec bhat = (ybar+tau);
  arma::vec logz = mat_digamma(ahat) - log(bhat);
  z = ahat/bhat;
  tau += sum(logz - z + log(tau) + 1 - R::digamma(tau))/(y.n_rows*(1/tau - R::trigamma(tau)));
}

double lowerbound_logML2(const arma::vec & z,
                         const arma::mat & alpha,
                         const arma::mat & beta,
                         const double & tau,
                         const arma::mat & lambda,
                         const arma::mat & loglambda,
                         const arma::vec & R,
                         const arma::vec & y,
                         const arma::uvec & xi,
                         const arma::uvec & xp,
                         const double & a,
                         const double & b){
  arma::vec ahat = (y+tau);
  arma::vec ybar = sum(myprod(y.n_rows, xi, xp, lambda), 1);
  arma::vec bhat = (ybar+tau);
  arma::vec logz = mat_digamma(ahat) - log(bhat);
  return sum(-R + y%log(R) - lgamma(y+1) + y%logz) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha))+
    + sum(z%bhat - z*a - ahat%log(bhat) + a*log(a)+lgamma(ahat)-lgamma(a));
}

// [[Rcpp::export]]
List doVB_negbin(const arma::vec & y,
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
  //arma::mat U = arma::zeros<arma::mat>(N,L);
  arma::vec z = arma::ones<arma::vec>(N);
  arma::vec ll(iter);
  double tau = 1.0;
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    up_B2(N, beta, V, logV, z, alpha, xi, xp, varind, b);
    up_tau(tau, z, y, V, xi,xp);
    ll.row(i) = lowerbound_logML2(z, alpha, beta, tau, V, logV, R, y, xi, xp, a, b);
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("precision")=tau,
                      Named("ELBO")=ll);
}
