#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_shape.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//rate parameters
void up_B(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
      arma::vec B = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k], varind[k+1]), vl) + b;
      beta.col(l).rows(varind[k], varind[k+1]-1) = B;
      V.col(l).rows(varind[k], varind[k+1]-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

double lowerbound_logML_pois(const arma::mat & alpha,
                        const arma::mat & beta,
                        const arma::mat & lambda,
                        const arma::mat & loglambda,
                        const arma::vec & R,
                        const arma::vec & y,
                        const double & a,
                        const double & b){
  return sum(-R + y%log(R) - lgamma(y+1)) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha));
}


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
