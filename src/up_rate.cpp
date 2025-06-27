#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate.h"
#include <random>
#include <chrono>
#include <iostream>
#include "TimerLogger.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

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
  for(int l = 0; l < L; l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl) + b;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
      V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

//geometric sampling for SVB
arma::vec geomsum(const int & D,
                  const arma::vec & Vl,
                  const arma::uvec & varind,
                  const double & rho,
                  const arma::uword & k,
                  const int & M_max){
  int n0 = R::rgeom(rho);
  arma::vec out = arma::zeros<arma::vec>(D);
  if (n0 == 0) {
    return out;
  }
  int M = std::min(M_max, n0);
  arma::uvec out_indices = arma::randi<arma::uvec>(M, arma::distr_param(0, D - 1));
  arma::vec fl(M);
  double MR = (double) n0 / (double) M;
  // sampling X with size n0
  for(int j = 0; j < M; j++){
    fl.fill(1.0);
    //product
    for(arma::uword i = 0; i < k; i++){
      arma::uvec indices;
      if (varind(i) == varind(i+1) - 1) {
        indices = arma::uvec(M, arma::fill::value(varind(i)));
      } else {
        const int start = varind(i);
        const int end = varind(i+1) - 1;
        indices = arma::randi<arma::uvec>(M, arma::distr_param(start, end));
      }
      fl %= Vl(indices);
    }
    //skip i==k
    for(arma::uword i = k + 1; i < (varind.n_rows - 1); i++){
      arma::uvec indices;
      if (varind(i) == varind(i+1) - 1) {
        indices = arma::uvec(M, arma::fill::value(varind(i)));
      } else {
        const int start = varind(i);
        const int end = varind(i+1) - 1;
        indices = arma::randi<arma::uvec>(M, arma::distr_param(start, end));
      }
      fl %= Vl(indices);
    }
    //sum i==k
    out(out_indices(j)) += MR * sum(fl);
  }
  return out;
}

// void up_Bs(const int & N,
//               arma::mat & beta,
//               arma::mat & V,
//               arma::mat & logV,
//               const arma::mat & alpha,
//               const arma::uvec & xi,
//               const arma::uvec & xp,
//               const arma::uvec & varind,
//               const double & N1,
//               const double & b,
//               const double & rho,
//               const int & M_max){
//   int K = varind.n_rows - 1;
//   int L = V.n_cols;
//   for(int l = 0; l < L; l++){
//     arma::vec vl = myprodvec(N, xi, xp, V.col(l));
//     for(int k=0; k < K; k++){
//       vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//       arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
//       arma::vec B0 = geomsum(varind(k+1) - varind(k), V.col(l), varind, rho, k, M_max);
//       arma::vec B = B1 + B0;
//       beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
//       V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
//       vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//     }
//     logV = mat_digamma(alpha) - log(beta);
//   }
// }

void up_Bs_sp(const int & N,
              arma::mat & beta,
              arma::mat & V,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const double & rho,
              const double & N1S,
              const double & NS,
              const double & b,
              const int & M_max){
  TimerLogger total_timer("upB_total");
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  arma::vec vl(N);
  for(int l = 0; l < L; l++){
    vl = myprodvec(N, xi, xp, V.col(l)); // # N
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l)); // # N
      arma::vec B1;
      //{
        //TimerLogger block1("B1");
        B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);        
      //}
      arma::vec B0;
      //{
        //TimerLogger block2("B0");
        B0 = geomsum(varind(k+1) - varind(k), V.col(l), varind, rho, k, M_max);
      //}
      arma::vec B = N1S*B1 + N1S*B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}

//initialize B
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

// void minus_Bs(const int & N,
//              arma::mat & beta,
//              arma::mat & V,
//              const arma::uvec & xi,
//              const arma::uvec & xp,
//              const arma::uvec & varind){
//   int K = varind.n_rows - 1;
//   int L = V.n_cols;
//   for(int l=0;l<L;l++){
//     arma::vec vl = myprodvec(N, xi, xp, V.col(l));
//     for(int k=0; k<K; k++){
//       vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//       arma::vec B = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
//       beta.col(l).rows(varind(k), varind(k+1) - 1) -= B;
//       vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//     }
//   }
// }