#include <RcppArmadillo.h>
#include "myproduct.h"
#include "logexpfuns.h"
#include "up_rate.h"
#include <random>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

int sample_index_cpp(int n) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, n - 1);
  return dis(gen);
}

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
                    const arma::mat & Vl,
                    const arma::uvec & varind,
                    const double & rho,
                    const arma::uword & k){
  arma::vec out = arma::zeros<arma::mat>(D);
  int M = R::rgeom(rho);
  // sampling X with size M
    for(int j = 0; j < M; j++){
      //product
      arma::vec fl = arma::ones<arma::vec>(M);
      for(arma::uword i = 0; i < k; i++){
        arma::distr_param support = arma::distr_param((int) varind(i), (int) varind(i+1) - 1);
        arma::uvec indices = arma::randi<arma::uvec>(M, support);
        fl %= Vl.rows(indices);
      }
      //skip i==k
      for(arma::uword i = k + 1; i < (varind.n_rows - 1); i++){
        arma::distr_param support = arma::distr_param((int) varind(i), (int) varind(i+1) - 1);
        arma::uvec indices = arma::randi<arma::uvec>(M, support);
        fl %= Vl.rows(indices);
      }
      //sum
      //i==k
      int index = sample_index_cpp(D);
      out.row(index) += sum(fl);
    }
  return out;
}

void up_Bs(const int & N,
              arma::mat & beta,
              arma::mat & V,
              arma::mat & logV,
              const arma::mat & alpha,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const double & N1,
              const double & b,
              const double & rho){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l = 0; l < L; l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      arma::vec B0 = geomsum(varind(k+1) - varind(k), V.col(l), varind, rho, k);
      arma::vec B = B1 + B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
      V.col(l).rows(varind(k), varind(k+1) - 1) = alpha.col(l).rows(varind(k), varind(k+1) - 1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

void up_Bs_sp(const int & N,
              arma::mat & beta,
              arma::mat & V,
              const arma::uvec & xi,
              const arma::uvec & xp,
              const arma::uvec & varind,
              const double & rho,
              const double & N1S, const double & NS,
              const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l = 0; l < L; l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k < K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
      arma::vec B1 = mysum_t(varind(k+1) - varind(k), xi, xp.rows(varind(k), varind(k+1)), vl);
      arma::vec B0 = geomsum(varind(k+1) - varind(k), V.col(l), varind, rho, k);
      arma::vec B = N1S*B1 + N1S*B0;
      beta.col(l).rows(varind(k), varind(k+1) - 1) = B + b;
      vl %= myprodvec_sub(N, xi, xp, varind(k), varind(k+1), V.col(l));
    }
  }
}
