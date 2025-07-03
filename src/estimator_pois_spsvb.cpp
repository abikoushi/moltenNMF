#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate.h"
#include "geometricsampling.h"
#include "ELBO.h"
#include "subset_sp.h"
#include "lr.h"
#include "rand.h"
#include <progress.hpp>
#include <progress_bar.hpp>
//#include "TimerLogger.h"
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

////
//update variational parameters
////
void up_theta_sp(const int & N,
                 arma::mat & alpha,
              arma::mat & beta,
              arma::vec & R,
              const arma::mat & V,
              const arma::mat & logV,
              const arma::vec & y,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const double & rho,
              const double & N1S,
              const double & a,
              const double & b,
              const int & M_max){
  //TimerLogger up_theta_total("up_theta");
  int K = (varind.n_rows - 1);
  int L = V.n_cols;
  arma::mat r =  myprod(y.n_rows, xi, xp, exp(logV));
  R = sum(r, 1);
  alpha = N1S*mysum_t(alpha.n_rows, xi, xp, r.each_col()%(y/R)) + a;
  arma::vec vl(y.n_rows);
  //for geometric sampling
  int n0 = R::rgeom(rho);
  int M = std::min(M_max, n0);
  double MR = (double) n0 / (double) M;
  arma::vec vl0(M);
  arma::umat U(M, K);
  ///
  for(int l = 0; l < L; l++){
    vl = r.col(l);
    vl0 = geomprod_all(M, V.col(l), varind, U);
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(y.n_rows, xi, xp, varind(k), varind(k+1), V.col(l)); // # N
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      //arma::vec B0 = geomsum(varind(k+1) - varind(k), V.col(l), varind, rho, k, M_max);
      arma::vec B0 = geomsum_k(varind(k+1) - varind(k), M, MR, vl0, V.col(l), varind, k, U);
      arma::vec B = N1S*B1 + N1S*B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
      vl %= myprodvec_sub(y.n_rows, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}

void upEV(arma::mat & V, arma::mat & logV,
          const arma::mat & alpha,
          const arma::mat & beta,
          const arma::uvec & up){
  //TimerLogger up_V_total("up_V");
  V.rows(up) = alpha.rows(up)/beta.rows(up);
  logV.rows(up) = mat_digamma(alpha.rows(up)) - log(beta.rows(up));
}

/////
//SVB main
/////
// [[Rcpp::export]]
List doSVB_pois_sp_skip(const int & N,
                    const arma::vec & yv,
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
                    const int & M_max,
                    const bool & display_progress){
  const int N1 =  yv.n_rows;
  const double N0 = N - N1;
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
      arma::uvec uid = bags.col(step);
      S_yv = yv.rows(uid);
      arma::uvec up;
      subset_spx(S_xi, S_xp, xi, xp, uid, up);
      arma::mat alpha_s = alpha;
      arma::mat beta_s = beta;
      arma::vec SR = R.rows(uid);
      up_theta_sp(N, alpha_s, beta_s,
                  SR, V, logV,
                  S_yv, S_xi, S_xp, varind,
                  p1, N1S, a, b, M_max);
      R(uid) = SR;
      alpha.rows(up) = rho2 * alpha.rows(up) + rho * alpha_s.rows(up);
      beta.rows(up) = rho2 * beta.rows(up) + rho * beta_s.rows(up);
      upEV(V, logV, alpha, beta, up);
    }
    ll.row(epoc) += lowerbound_logML_pois(alpha, beta, V, logV, R, yv, a, b);      
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}