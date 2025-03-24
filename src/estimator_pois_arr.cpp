#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include "logexpfuns.h"
#include "rand.h"
#include "outerprod.h"
#include "up_shape_2D.h"
#include "up_param_array.h"
#include "lr.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

//without weight
double up_B_arr(const arma::field<arma::mat> & alpha,
            arma::mat & beta,
            arma::field<arma::mat> & V,
            arma::field<arma::mat> & logV,
            const double & b,
            const int & L,
            const arma::uvec & dims, 
            const int & k){
  double lp = 0;
  //const int K = V.n_rows;
  for(int l=0; l<L; l++){
    double B0 = sumouterprod(V, k, l);
    lp -= B0;
    beta(k,l) = B0 + b;
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}

//with weight
double up_B_arr(const arma::field<arma::mat> & alpha,
            arma::mat & beta,
            arma::field<arma::mat> & V,
            arma::field<arma::mat> & logV,
            const double & b,
            const int & L,
            const arma::uvec & dims, 
            const arma::field<arma::vec> & weight,
            const int & k){
  double lp = 0;
  //const int K = V.n_rows;
  for(int l=0; l<L; l++){
    double B0 = sumouterprod(V, k, l, weight);
    lp -= B0;
    beta(k,l) = B0 + b;
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}

//without weight
double up_theta_arr(arma::field<arma::mat> & alpha,
                arma::mat & beta,
                arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const int & L,
                const arma::uvec & dims, 
                const double & a,
                const double & b){
  double lp_a = 0;
  double lp_b = 0;
  for(int k=0; k<dims.n_rows; k++){
    lp_a += up_A_2D(alpha, logV, y, X, a, L, dims, k);
    lp_b += up_B_arr(alpha, beta, V, logV, b, L, dims,k);    
  }
  return lp_a+lp_b;
}

//with weight
double up_theta_arr(arma::field<arma::mat> & alpha,
                arma::mat & beta,
                arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const int & L,
                const arma::uvec & dims, 
                const double & a,
                const double & b,
                const arma::field<arma::vec> & weight){
  double lp_a = 0;
  double lp_b = 0;
  for(int k=0; k<dims.n_rows; k++){
    lp_a += up_A_2D(alpha, logV, y, X, a, L, dims, k);
    lp_b += up_B_arr(alpha, beta, V, logV, b, L, dims, weight, k);
  }
  return lp_a+lp_b;
}

// [[Rcpp::export]]
List doVB_pois_arr(const arma::vec & y,
               const arma::umat & X,
               const arma::uvec & dims,
               const int & L,
               const int & iter,
               const double & a,
               const double & b, 
               const bool & display_progress){
  arma::field<arma::mat> V(dims.n_rows);
  arma::field<arma::mat> logV(dims.n_rows);
  arma::field<arma::mat> alpha(dims.n_rows);
  arma::mat beta(dims.n_rows, L);
  for(int k=0; k<dims.n_rows; k++){
    V(k) = randg<mat>(dims(k),L);
    logV(k) = log(V(k));
    beta.row(k) = sum(V(k), 0);
  }
  beta += b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    double lp0 = up_theta_arr(alpha, beta, V, logV, y, X, L, dims, a, b);
    lp(i) = lp0 + kld2(alpha, beta, a, b);
    pb.increment();
  }
  lp -= sum(lgamma(y+1)); 
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("logprob")=lp);
}

// [[Rcpp::export]]
List doVB_pois_w_arr(const arma::vec & y,
                 const arma::umat & X,
                 const arma::uvec & dims,
                 const int & L,
                 const int & iter,
                 const double & a,
                 const double & b,
                 const arma::field<arma::vec> & weight,
                 const bool & display_progress){
  field<mat> V(X.n_cols);
  field<mat> logV(X.n_cols);
  field<mat> alpha(X.n_cols);
  arma::mat beta(X.n_cols, L);
  for(int k=0; k<X.n_cols; k++){
    V(k) = randg<mat>(dims(k),L);
    logV(k) = log(V(k));
    beta.row(k) = sum(V(k), 0);
  }
  beta += b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    double lp0 = up_theta_arr(alpha, beta, V, logV, y, X, L, dims, a, b, weight);
    lp(i) = lp0 + kld2(alpha, beta, a, b);
    pb.increment();
  }
  lp -= sum(lgamma(y+1)); 
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("logprob")=lp);
}


////
//svb
//stochastic mini-batches
//

//without weight & use SumV
double up_Bs_arr(const arma::field<arma::mat> & alpha,
             arma::mat & beta,
             arma::field<arma::mat> & V,
             arma::field<arma::mat> & logV,
             arma::rowvec & sumVk,
             const arma::field<arma::uvec> & uid,
             const double & b,
             const int & L,
             const double & NS,
             const int & k){
  double lp = 0;
  for(int l=0; l<L; l++){
    sumVk(l) += NS*sumouterprod_s(V, k, l, uid);
    lp -= sumVk(l) - b;
    beta(k,l) = sumVk(l);
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}

double up_Bs_arr(const arma::field<arma::mat> & alpha,
                 arma::mat & beta,
                 arma::field<arma::mat> & V,
                 arma::field<arma::mat> & logV,
                 arma::rowvec & sumVk,
                 const arma::field<arma::uvec> & uid,
                 const double & b,
                 const int & L,
                 const int & k){
  double lp = 0;
  for(int l=0; l<L; l++){
    sumVk(l) += sumouterprod_s(V, k, l, uid);
    lp -= sumVk(l) - b;
    beta(k,l) = sumVk(l);
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}

//without weight & use SumV
double up_theta_s_arr(arma::field<arma::mat> & alpha,
                  arma::mat & beta,
                  arma::field<arma::mat> & V,
                  arma::field<arma::mat> & logV,
                  const arma::vec & y,
                  const arma::umat & X,
                  const int & L,
                  const arma::field<arma::uvec> & uid,
                  const double & a,
                  const double & b,
                  const double & NS){
  double lp_a = 0;
  double lp_b = 0;
  for(int k=0; k<X.n_cols; k++){
    //Rprintf("%d : ", k);
    lp_a += up_As_2D(alpha, logV, y, X, a, L, uid, k, NS);
    arma::rowvec sumVk = beta.row(k);
    for(int l = 0; l < L; l++){
      //sumVk(l) -= NS*sumouterprod_s(V, k, l, uid);
      sumVk(l) -= sumouterprod_s(V, k, l, uid);
    }
    lp_b += up_Bs_arr(alpha, beta, V, logV, sumVk, uid, b, L, k);
  }
  //Rprintf("\n");
  return lp_a+lp_b;
}


double doVB_pois_s_sub_arr(const arma::vec & y,
                       const arma::umat & X,
                       const arma::uvec & dims,
                       const int & L,
                       const int & iter,
                       const double & a,
                       const double & b,
                       const double & N1, 
                       const double & NS, 
                       const arma::field<arma::uvec> uid,
                       arma::field<arma::mat> & alpha, 
                       arma::mat & beta,
                       arma::field<arma::mat> & V,
                       arma::field<arma::mat> & logV){
  double lp = 0.0;
  for (int i=0; i<iter; i++) {
    double lp0 = up_theta_s_arr(alpha, beta, V, logV, y, X, L, uid, a, b, NS);
    lp = lp0 + kld2(alpha, beta, a, b);
  }
  return lp;
}

//SVB without weight
// [[Rcpp::export]]
List doVB_pois_s_arr(const arma::umat & X,
                 const arma::vec & y,
                 const arma::uvec & dims,
                 const int & L,
                 const int & iter,
                 const int & subiter,
                 const double & a, const double & b,
                 const double & N1,
                 const int & bsize,
                 const arma::vec & lr_param,
                 const std::string & lr_type,
                 const bool & display_progress){
  field<mat> V(dims.n_rows);
  field<mat> logV(dims.n_rows);
  field<mat> alpha(dims.n_rows);
  for(int k=0; k<dims.n_rows; k++){
    alpha(k) = arma::ones<arma::mat>(dims(k), L);
  }
  arma::mat beta(dims.n_rows, L);
  for(int k=0; k<dims.n_rows; k++){
    arma::mat Vk = arma::mat(dims(k),L);
    V(k) = arma::randg<mat>(dims(k), L);
    logV(k) = log(V(k));
    beta.row(k) = sum(V(k), 0);
  }
  beta += b;
  std::unique_ptr<lr> g;
  if(lr_type == "const"){
    g.reset(new lr_const);
  }else if(lr_type == "exponential"){
    g.reset(new lr_power);
  }else{
    Rcpp::stop("This lr_type is not implemented\n");
  }
  // const double NS = (bsize/N1);
  const double NS = N1 / ((double) bsize);
  double invS = 1.0 / ( (double) bsize ); 
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::vec val(bsize);
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    rho *= invS;
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      arma::umat SX = X.rows(bag);
      arma::vec Sy = y.rows(bag);
      arma::field<arma::uvec> uid(dims.n_rows);
      for(int j = 0; j < X.n_cols; j++){
        uid(j) = unique(SX.col(j));
      }
      arma::field<arma::mat> alpha_s = alpha;
      arma::mat beta_s = beta;
      lp(epoc) += doVB_pois_s_sub_arr(Sy, SX, dims, L, subiter, a, b, N1, NS, uid, alpha_s, beta_s, V, logV);
      up_vpar(rho, rho2, uid, alpha, alpha_s, beta, beta_s);
    }
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("logprob")=lp);
}

//ToDo : 
//SVB with weight
//with weight & use sumV
double up_Bs_arr(const arma::field<arma::mat> & alpha,
                 arma::mat & beta,
                 arma::field<arma::mat> & V,
                 arma::field<arma::mat> & logV,
                 arma::rowvec & sumVk,
                 const arma::field<arma::uvec> & uid,
                 const double & b,
                 const int & L,
                 const double & NS,
                 const arma::field<arma::vec> & weight,
                 const int & k){
  double lp = 0;
  for(int l=0; l<L; l++){
    sumVk(l) += NS*sumouterprod_s(V, k, l, weight, uid);
    lp -= sumVk(l) - b;
    beta(k,l) = sumVk(l);
    arma::vec alpha_k_l = alpha(k).col(l);
    arma::vec V_l = alpha_k_l/beta(k,l);
    arma::mat Vk = V(k);
    Vk.col(l) = V_l;
    V(k) = Vk;
    arma::mat logv = logV(k);
    up_log_gamma(logv, alpha_k_l, log(beta(k,l)), l);
    logV(k) = logv;
  }
  return lp;
}

//with weight
double up_theta_s_arr(arma::field<arma::mat> & alpha,
                      arma::mat & beta,
                      arma::field<arma::mat> & V,
                      arma::field<arma::mat> & logV,
                      const arma::vec & y,
                      const arma::umat & X,
                      const int & L,
                      const arma::field<arma::uvec> & uid,
                      const double & a,
                      const double & b,
                      const double & NS,
                      const arma::field<arma::vec> & weight){
  double lp_a = 0;
  double lp_b = 0;
  for(int k=0; k<X.n_cols; k++){
    lp_a += up_As_2D(alpha, logV, y, X, a, L, uid, k);
    arma::rowvec sumVk = beta.row(k) - b;
    for(int l=0; l<sumVk.n_elem; l++){
      sumVk(l) -=  NS*sumouterprod_s(V, k, l, uid);      
    }
    lp_b += up_Bs_arr(alpha, beta, V, logV, sumVk, uid, b, L, NS, weight, k);
  }
  return lp_a+lp_b;
}