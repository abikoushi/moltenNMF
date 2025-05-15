#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate.h"
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

arma::vec ZeroSampling(const arma::vec & probX0,
                       const arma::uvec & varind,
                       const arma::vec & vl,
                       int k){
  arma::vec b0 = probX0.rows(varind(k), varind(k+1)-1)%vl.rows(varind(k), varind(k+1)-1);
  return b0;
}

void up_B_sp(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const arma::vec & probX0,
          const double & N0,
          const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l = 0; l < L; l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl) + b;
      arma::vec B0 = ZeroSampling(probX0, varind, vl, k);
      B0 *= N0;
      arma::vec B = B1 + B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
      V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

arma::vec randomsum(const int & D,
                    const arma::mat & Vl,
                    const arma::uvec & varind,
                    const arma::uvec & xp0,
                    const int & k){
  int count = xp0(k+1)-xp0(k);
  arma::vec fl = arma::ones<arma::vec>(count);
  for(arma::uword i = 0; i < varind.n_rows - 1; i++){
    arma::uvec indices = arma::randi<arma::uvec>(count, arma::distr_param((int) varind(i), (int) varind(i+1) - 1));
    fl %= Vl(indices);
  }
  return mysum_t_rv(D, xp0.rows(varind(k), varind(k+1)), fl);
}


void up_B_sp2(const int & N,
             arma::mat & beta,
             arma::mat & V,
             arma::mat & logV,
             const arma::mat & alpha,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const arma::uvec & varind,
             const arma::uvec & xp0,
             const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l = 0; l < L; l++){
    arma::vec fl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k < K; k++){
      fl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), fl) + b;
      //arma::vec B0 = mysum_t_rv(varind(k+1) - varind(k), xp0.rows(varind(k), varind(k+1)), fl);
      arma::vec B0 = randomsum(varind(k+1) - varind(k), V.col(l), varind, xp0, k);
      arma::vec B = B1 + B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
      V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
      fl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

void up_Bs(const int & N,
           arma::mat & beta,
           arma::mat & V,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl) + b;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}


void up_Bs_sp(const int & N,
             arma::mat & beta,
             arma::mat & V,
             const arma::uvec & xi,
             const arma::uvec & xp,
             const arma::uvec & varind,
             const arma::vec & probX0,
             const double & N0, const double & NS,
             const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l = 0; l < L; l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      arma::vec B0 = ZeroSampling(probX0, varind, vl, k);
      arma::vec B = NS*B1 + N0*B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}

// void up_B_sp(const int & N,
//           arma::mat & beta,
//           arma::mat & V,
//           arma::mat & logV,
//           const arma::mat & alpha,
//           const arma::uvec & xi,
//           const arma::uvec & xp,
//           const arma::uvec & varind,
//           const arma::vec & weight,
//           const double & b){
//   int K = varind.n_rows - 1;
//   int L = V.n_cols;
//   for(int l = 0; l < L; l++){
//     arma::vec vl = myprodvec(N, xi, xp, V.col(l));
//     for(int k=0; k<K; k++){
//       vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//       arma::vec B = mysum_t_ww(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl, weight) + b;
//       beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
//       V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
//       vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//     }
//     logV = mat_digamma(alpha) - log(beta);
//   }
// }

// void up_B_ww(const int & N,
//           arma::mat & beta,
//           arma::mat & V,
//           arma::mat & logV,
//           const arma::mat & alpha,
//           const arma::uvec & xi,
//           const arma::uvec & xp,
//           const arma::uvec & varind,
//           const arma::vec & weight,
//           const double & b){
//   int K = varind.n_rows - 1;
//   int L = V.n_cols;
//   for(int l = 0; l < L; l++){
//     arma::vec vl = myprodvec(N, xi, xp, V.col(l));
//     for(int k=0; k<K; k++){
//       vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//       arma::vec B = mysum_t_ww(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl, weight) + b;
//       beta.col(l).rows(varind(k), varind(k+1) - 1) = B;
//       V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
//       vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
//     }
//     logV = mat_digamma(alpha) - log(beta);
//   }
// }
