#include "RcppArmadillo.h"
#include "up_param_array.h"
#include "logexpfuns.h"
// [[Rcpp::depends(RcppArmadillo)]]

//update variational posterior
void up_vpar(const double rho, const double rho2,
             const arma::field<arma::uvec> uid,
             arma::field<arma::mat> & alpha,
             arma::field<arma::mat> & alpha_s,
             arma::mat & beta, arma::mat & beta_s){
  int K = beta.n_rows; 
  for(int j = 0; j < K; j++){
    alpha(j).rows(uid(j)) = rho2 * alpha(j).rows(uid(j)) + rho * alpha_s(j).rows(uid(j));
  }
  beta = rho2 * beta + rho * beta_s;  
}

//update latent V & logV
void up_latentV(arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::mat> & alpha,
                const arma::mat & beta, 
                const int & k){
  arma::mat Vk = V(k);
  arma::mat logv = logV(k);
  arma::mat alpha_k = alpha(k);
  for(int l=0; l < beta.n_cols; l++){
    arma::vec alpha_k_l = alpha_k.col(l);
    Vk.col(l) = alpha_k_l/beta(k,l);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
  }
  logV(k) = logv;
  V(k) = Vk;
}