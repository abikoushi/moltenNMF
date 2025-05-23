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

/////
//SVB
/////
//rate parameters
void plus_Bs(const int & N,
          arma::mat & beta,
          arma::mat & V,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      beta.col(l).rows(varind(k), varind(k+1) - 1) += B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}

void minus_Bs(const int & N,
             arma::mat & beta,
             arma::mat & V,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const arma::uvec & varind){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      beta.col(l).rows(varind(k), varind(k+1) - 1) -= B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}

void upEV(arma::mat & V, arma::mat & logV,
          const arma::mat & alpha, const arma::mat & beta){
  V = alpha/beta;
  logV = mat_digamma(alpha) - log(beta);
}

void upR(arma::vec & R,
         const arma::mat & loglambda,
         const arma::uvec & xi,
         const arma::uvec & xp,
         const arma::uvec uid){
    arma::mat r =  myprod(R.n_rows, xi, xp, exp(loglambda)); //(N, L)
    R.rows(uid) = sum(r.rows(uid), 1);
} 


// [[Rcpp::export]]
List doSVB_pois_sp_skip(const int & N,
                    const arma::vec & yv,
                    const arma::uvec & xi,
                    const arma::uvec & xp,
                    const arma::uvec & varind,
                    const double & N0,
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
  int N1 =  yv.n_rows;
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta(D, L);
  beta.fill(b);
  plus_Bs(N, beta, V, xi, xp, varind);
  arma::vec R = arma::zeros<arma::vec>(N1);
  arma::vec ll = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  const double N1S = ((double) N1) / ((double) bsize);
  const double NS = ((double) N) / ((double) bsize);
  const double p1 = ((double) N1) / ((double) N);
  Progress pb(iter, display_progress);
  for(int epoc = 0; epoc < iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    for(int step = 0; step < (int) bags.n_cols; step++){
      arma::vec S_yv;
      arma::uvec S_xp;
      arma::uvec S_xi;
      S_yv = yv.rows(bags.col(step));
      subset_spx(S_xi, S_xp, xi, xp, bags.col(step));
      arma::mat alpha_s = alpha;
      arma::mat beta_s = beta;
      arma::vec SR = R.rows(bags.col(step));
      up_As_sp2(alpha_s, SR, logV, S_yv, S_xi, S_xp, a, N1S);
      up_Bs_sp(N, beta_s, V, S_xi, S_xp, varind, p1, N1S, NS, b);
      upEV(V, logV, alpha, beta);
      alpha = rho2 * alpha + rho * alpha_s;
      beta = rho2 * beta + rho * beta_s;
      upR(R, logV, xi, xp, bags.col(step));
    }
    ll.row(epoc) += lowerbound_logML_pois(alpha, beta, V, logV, R, yv, a, b);      
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}