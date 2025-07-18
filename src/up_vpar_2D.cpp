#include "RcppArmadillo.h"
#include "up_vpar_2D.h"
#include "logexpfuns.h"
// [[Rcpp::depends(RcppArmadillo)]]

//update V & logV
void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta){
  for(int k = 0; k < (int) beta.n_rows; k++){
    arma::mat Vk = V(k);
    arma::mat logv = logV(k);
    arma::mat alpha_k = alpha(k);
    for(int l = 0; l < (int) beta.n_cols; l++){
      arma::vec alpha_k_l = alpha_k.col(l);
      Vk.col(l) = alpha_k_l/beta(k,l);
      up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    }
    logV(k) = logv;
    V(k) = Vk;
  }
}

void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta, 
                const int & k){
  arma::mat Vk = V(k);
  arma::mat logv = logV(k);
  arma::mat alpha_k = alpha(k);
  for(int l = 0; l < (int) beta.n_cols; l++){
    arma::vec alpha_k_l = alpha_k.col(l);
    Vk.col(l) = alpha_k_l/beta(k,l);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
  }
  logV(k) = logv;
  V(k) = Vk;
}

void up_latentV_uid(arma::field<arma::mat> & V,
                    arma::field<arma::mat> & logV,
                    const arma::field<arma::mat> & alpha,
                    const arma::mat & beta, 
                    const int & k,
                    const arma::field<arma::uvec> & uid){
  arma::mat Vk = V(k);
  arma::mat logv = logV(k);
  arma::mat alpha_k = alpha(k);
  arma::uvec uid_k = uid(k);
  for(int l = 0; l < (int) beta.n_cols; l++){
    arma::vec alpha_k_l = alpha_k.col(l);
    arma::vec vl = Vk.col(l);
    vl.rows(uid_k) = alpha_k_l.rows(uid_k)/beta(k,l);
    Vk.col(l) = vl;
    up_log_gamma_uid(logv, alpha_k_l, log(beta(k,l)), l, uid_k);
  }
  logV(k) = logv;
  V(k) = Vk;
}


arma::rowvec sumV_uid(const arma::field<arma::mat> V, 
                      const arma::field<arma::uvec> uid,
                      const int k){
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  arma::rowvec SumV = sum(V(not_k).rows(uid(not_k)), 0);
  return SumV;
}

arma::rowvec sumV_uid(const arma::field<arma::mat> V, 
                      const arma::field<arma::uvec> uid,
                      const int k,
                      const arma::field<arma::vec> weight){
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  arma::mat Vk = V(not_k).rows(uid(not_k));
  Vk.each_col() %= weight(k).rows(uid(not_k));
  arma::rowvec SumV = sum(Vk, 0);
  return SumV;
}

void minus_sumV(arma::mat & beta, 
                const arma::field<arma::mat> V, 
                const arma::field<arma::uvec> uid,
                const int k){
  beta.row(k) -= sumV_uid(V, uid, k);
}

void plus_sumV(arma::mat & beta, 
               const arma::field<arma::mat> V, 
               const arma::field<arma::uvec> uid,
               const int k){
  beta.row(k) += sumV_uid(V, uid, k);
}

void minus_sumV(arma::mat & beta, 
                const arma::field<arma::mat> V, 
                const arma::field<arma::uvec> uid,
                const int k,
                const arma::field<arma::vec> weight){
  beta.row(k) -= sumV_uid(V, uid, k, weight);
}

//update variational posterior
void up_vpar(arma::field<arma::mat> & alpha,
             arma::mat & beta,
             arma::field<arma::mat> & V,
             arma::field<arma::mat> & logV,
             const arma::field<arma::mat> & alpha_s,
             const arma::mat & beta_s,
             const arma::field<arma::uvec> & uid,
             const double rho,
             const double & rho2){
  int K = beta.n_rows;
  beta = rho2 * beta + rho * beta_s;
  for(int j = 0; j < K; j++){
    alpha(j).rows(uid(j)) = rho2 * alpha(j).rows(uid(j)) + rho * alpha_s(j).rows(uid(j));
    up_latentV_uid(V, logV, alpha, beta, j, uid);
  }
}
