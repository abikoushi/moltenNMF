#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "logexpfuns.h"
#include "ELBO.h"
#include "rand.h"
#include "up_shape_2D.h"
#include "lr.h"
#include "up_param_array.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

arma::rowvec sumV(const arma::field<arma::mat> V, 
                  const int k){
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  arma::rowvec SumV = sum(V(not_k), 0);
  return SumV;
}

arma::rowvec sumV(const arma::field<arma::mat> V, 
                  const int k,
                  const arma::field<arma::vec> weight){
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  arma::mat Vk = V(not_k);
  Vk.each_col() %= weight(not_k);
  arma::rowvec SumV = sum(Vk, 0);
  return SumV;
}

double up_B_2D(const arma::field<arma::mat> & alpha,
               arma::mat & beta,
               arma::field<arma::mat> & V,
               arma::field<arma::mat> & logV,
               const double & b,
               const int & L,
               const arma::uvec & dims, 
               const int & k){
  double lp = 0;
  arma::rowvec B0 = sumV(V, k);
  lp -= sum(B0);
  beta.row(k) = B0 + b;
  up_latentV(V, logV, alpha, beta, k);
  return lp;
}

void up_B_2D(arma::mat & beta,
               arma::field<arma::mat> & V,
               const double & b){
  //double lp = 0;
  //row k = 0; column k = 1
  for(int k=0; k<beta.n_rows ;k++){
    arma::rowvec B0 = sumV(V, k);
    beta.row(k) = B0 + b;    
    //lp -= sum(B0);
  }
}

double up_theta_2D(arma::field<arma::mat> & alpha,
                   arma::mat & beta,
                   arma::field<arma::mat> & V,
                   arma::field<arma::mat> & logV,
                   const arma::vec & y,
                   const arma::umat & X,
                   const int & L,
                   const arma::uvec & dims, 
                   const double & a,
                   const double & b){
  double lp = 0;
  //lp += up_A_2D(alpha, logV, y, X, a, L, dims);
  //up_B_2D(beta, V, b);
  //row k = 0
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 0);
  up_B_2D(alpha, beta, V, logV, b, L, dims, 0);
  //column k = 1
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 1);
  up_B_2D(alpha, beta, V, logV, b, L,dims, 1);
  return lp;
}

// [[Rcpp::export]]
List doVB_pois_2D(arma::field<arma::mat> V,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::uvec & dims,
                  const int & L,
                  const int & iter,
                  const double & a,
                  const double & b, 
                  const bool & display_progress){
  arma::umat X(y.n_rows, 2);
  X.col(0) = rowi;
  X.col(1) = coli;
  //arma::field<arma::mat> V(2);
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  V(0) = arma::randg<arma::mat>(dims(0),L);
  V(1) = arma::randg<arma::mat>(dims(1),L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(dims(0), L);
  alpha(1) = arma::ones<arma::mat>(dims(1), L);
  arma::mat beta = arma::ones<arma::mat>(2, L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    lp(i) = up_theta_2D(alpha, beta, V, logV, y, X, L, dims, a, b);
    lp(i) += kld2(alpha, beta, a, b);
    up_latentV(V, logV, alpha, beta);
    pb.increment();
  }
  lp -= sum(lgamma(y+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}

//////
//with weight, batch
//////
double up_B_2D(const arma::field<arma::mat> & alpha,
               arma::mat & beta,
               arma::field<arma::mat> & V,
               arma::field<arma::mat> & logV,
               const arma::field<arma::vec> & weight,
               const double & b,
               const int & L,
               const arma::uvec & dims, 
               const int & k){
  double lp = 0;
  arma::rowvec B0 = sumV(V, k, weight);
  lp -= sum(B0);
  beta.row(k) = B0 + b;
  up_latentV(V, logV, alpha, beta, k);
  return lp;
}


double up_theta_2D(arma::field<arma::mat> & alpha,
                   arma::mat & beta,
                   arma::field<arma::mat> & V,
                   arma::field<arma::mat> & logV,
                   const arma::field<arma::vec> & weight,
                   const arma::vec & y,
                   const arma::umat & X,
                   const int & L,
                   const arma::uvec & dims, 
                   const double & a,
                   const double & b){
  double lp = 0;
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 0);
  lp += up_B_2D(alpha, beta, V, logV, weight, b, L, dims, 0);
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 1);
  lp += up_B_2D(alpha, beta, V, logV, weight, b, L,dims, 1);
  return lp;
}

// [[Rcpp::export]]
List doVB_pois_2D_ww(arma::field<arma::mat> V,
                  const arma::vec & y,
                  const arma::uvec & rowi,
                  const arma::uvec & coli,
                  const arma::uvec & dims,
                  const int & L,
                  const int & iter,
                  const arma::field<arma::vec> & weight,
                  const double & a,
                  const double & b, 
                  const bool & display_progress){
  arma::umat X(y.n_rows, 2);
  X.col(0) = rowi;
  X.col(1) = coli;
  //arma::field<arma::mat> V(2);
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  V(0) = arma::randg<arma::mat>(dims(0),L);
  V(1) = arma::randg<arma::mat>(dims(1),L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(dims(0), L);
  alpha(1) = arma::ones<arma::mat>(dims(1), L);
  arma::mat beta = arma::ones<arma::mat>(2, L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    double lp0 = up_theta_2D(alpha, beta, V, logV, weight, y, X, L, dims, a, b);
    lp(i) = lp0 + kld2(alpha, beta, a, b);
    pb.increment();
  }
  lp -= sum(lgamma(y+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}


//////
//SVB, without weight
//////
//without weight
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

void minus_sumV(arma::mat & beta, 
                const arma::field<arma::mat> V, 
                const arma::field<arma::uvec> uid,
                const int k,
                const double b){
  beta.row(k) -= sumV_uid(V, uid, k);
}

void plus_sumV(arma::mat & beta, 
               const arma::field<arma::mat> V, 
               const arma::field<arma::uvec> uid,
               const int k,
               const double b){
  beta.row(k) += sumV_uid(V, uid, k);
}


double up_Bs_2D(const arma::field<arma::mat> & alpha,
                arma::mat & beta,
                arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::uvec> & uid,
                const double & b,
                const int & L,
                const int & k){
  double lp = 0;
  beta.row(k) = sumV(V, k) + b;
  //plus_sumV(beta, V, uid, k, b);
  up_latentV_uid(V, logV, alpha, beta, k, uid);
  lp -= sum((beta.row(k) - b));
  return lp;
}

double up_theta_s_2D(arma::field<arma::mat> & alpha,
                     arma::mat & beta,
                     arma::field<arma::mat> & V,
                     arma::field<arma::mat> & logV,
                     const int & L,
                     const arma::vec & y,
                     const arma::umat & X,
                     const arma::field<arma::uvec> & uid,
                     const double & a,
                     const double & b,
                     const double & NS){
  double lp = 0;
  //rows k = 0
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 0, NS);
  //minus_sumV(beta, V, uid, 0, b);
  lp += up_Bs_2D(alpha, beta, V, logV, uid, b, L, 0);
  //columns k = 1
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 1, NS);
  //minus_sumV(beta, V, uid, 1, b);
  lp += up_Bs_2D(alpha, beta, V, logV, uid, b, L,  1);
  return lp;
}

double doVB_pois_s_sub_2D(const arma::vec & y,
                          const arma::umat & X,
                          const int & L,
                          const int & iter,
                          const double & a,
                          const double & b,
                          const double & NS, 
                          const arma::field<arma::uvec> & uid,
                          arma::field<arma::mat> & alpha, 
                          arma::mat & beta,
                          arma::field<arma::mat> & V,
                          arma::field<arma::mat> & logV){
  double lp = 0.0;
  for (int i=0; i<iter; i++) {
    double lp1 = up_theta_s_2D(alpha, beta, V, logV, L, y, X, uid, a, b, NS);
    lp = lp1 + kld2(alpha, beta, a, b);
  }
  return lp;
}

//without weight
// [[Rcpp::export]]
List doVB_pois_s_2D(const arma::vec & y,
                    const arma::uvec & rowi,
                    const arma::uvec & coli,
                    const int & L,
                    const int & iter,
                    const int & subiter,
                    const double & a,
                    const double & b,
                    const double & N1,
                    const int & Nr,
                    const int & Nc,
                    const int & bsize,
                    const arma::vec & lr_param,
                    const std::string & lr_type,
                    const bool & display_progress){
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  arma::field<arma::mat> V(2);
  V(0) = arma::randg<arma::mat>(Nr,L);
  V(1) = arma::randg<arma::mat>(Nc,L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(Nr, L);
  alpha(1) = arma::ones<arma::mat>(Nc, L);
  arma::mat beta(2,L);
  beta.row(0) = sum(V(1), 0) + b;
  beta.row(1) = sum(V(0), 0) + b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  const double NS = N1 / ((double) bsize);
  double invS = 1.0 / ( (double) bsize ); 
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    //rho *= invS;
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      arma::vec Sy = y.rows(bag);
      arma::umat SX(bag.n_rows, 2);
      SX.col(0) = rowi.rows(bag);
      SX.col(1) = coli.rows(bag);
      arma::field<arma::uvec> uid(2);
      uid(0) = unique(SX.col(0));
      uid(1) = unique(SX.col(1));
      arma::field<arma::mat> alpha_s = alpha;
      arma::mat beta_s = beta;
      lp(epoc) += up_As_2D(alpha_s, logV, Sy, SX, a, L, uid, NS);
      up_B_2D(beta_s,V,b);
      //lp(epoc) += doVB_pois_s_sub_2D(Sy, SX, L, subiter, a, b, NS, uid, alpha_s, beta_s, V, logV);
      up_vpar(rho, rho2, uid, alpha, alpha_s, beta, beta_s, V, logV);
    }
    lp(epoc) += kld2(alpha, beta, a, b);
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}

// use only 1 sample in each updates
// [[Rcpp::export]]
List doVB_pois_s_2D_t1(const arma::vec & y,
                    const arma::uvec & rowi,
                    const arma::uvec & coli,
                    const int & L,
                    const int & iter,
                    const int & subiter,
                    const double & a,
                    const double & b,
                    const double & N1,
                    const int & Nr,
                    const int & Nc,
                    const int & bsize,
                    const arma::vec & lr_param,
                    const std::string & lr_type,
                    const bool & display_progress){
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  arma::field<arma::mat> V(2);
  V(0) = arma::randg<arma::mat>(Nr,L);
  V(1) = arma::randg<arma::mat>(Nc,L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(Nr, L);
  alpha(1) = arma::ones<arma::mat>(Nc, L);
  arma::mat beta(2,L);
  beta.row(0) = sum(V(1), 0) + b;
  beta.row(1) = sum(V(0), 0) + b;
  arma::vec lp = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  const double NS = N1 / ((double) bsize);
  // double invS = 1.0 / ( (double) bsize ); 
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    // rho *= invS;
    arma::uvec bag = arma::randperm(N1, bsize);
      arma::vec Sy = y.rows(bag);
      arma::umat SX(bag.n_rows, 2);
      SX.col(0) = rowi.rows(bag);
      SX.col(1) = coli.rows(bag);
      arma::field<arma::uvec> uid(2);
      uid(0) = unique(SX.col(0));
      uid(1) = unique(SX.col(1));
      arma::field<arma::mat> alpha_s = alpha;
      arma::mat beta_s = beta;
      lp(epoc) += doVB_pois_s_sub_2D(Sy, SX, L, subiter, a, b, NS, uid, alpha_s, beta_s, V, logV);
      up_vpar(rho, rho2, uid, alpha, alpha_s, beta, beta_s);
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}

////////
///with weight, SVB
////////
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
                const int k,
                const double b,
                const arma::field<arma::vec> weight){
  beta.row(k) -= sumV_uid(V, uid, k, weight);
  //beta.elem(find(beta<0.0)).fill(b);
}

double up_Bs_2D(const arma::field<arma::mat> & alpha,
                arma::mat & beta,
                arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                const arma::field<arma::vec> & weight,
                const arma::field<arma::uvec> & uid,
                const double & b,
                const int & L,
                const int & k){
  double lp = 0;
  beta.row(k) = sumV_uid(V, uid, k, weight);
  lp -= sum((beta.row(k) - b));
  up_latentV_uid(V, logV, alpha, beta, k, uid);
  return lp;
}

double up_theta_s_2D(arma::field<arma::mat> & alpha,
                     arma::mat & beta,
                     arma::field<arma::mat> & V,
                     arma::field<arma::mat> & logV,
                     const arma::field<arma::vec> & weight,
                     const int & L,
                     const arma::vec & y,
                     const arma::umat & X,
                     const arma::field<arma::uvec> & uid,
                     const double & a,
                     const double & b,
                     const double & NS){
  double lp = 0;
  //rows k=0
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 0, NS);
  lp += up_Bs_2D(alpha, beta, V, logV, weight, uid, b, L, 0);
  //columns k=1
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 1, NS);
  lp += up_Bs_2D(alpha, beta, V, logV, weight, uid, b, L, 1);
  return lp;
}

double doVB_pois_s_sub_2D(const arma::vec & y,
                          const arma::umat & X,
                          const int & L,
                          const int & iter,
                          const double & a,
                          const double & b,
                          const double & NS, 
                          const arma::field<arma::uvec> & uid,
                          const arma::field<arma::vec> & weight,
                          arma::field<arma::mat> & alpha, 
                          arma::mat & beta,
                          arma::field<arma::mat> & V,
                          arma::field<arma::mat> & logV){
  double lp = 0.0;
  for (int i=0; i<iter; i++) {
    double lp1 = up_theta_s_2D(alpha, beta, V, logV, weight, L, y, X, uid, a, b, NS);
    lp = lp1 + kld2(alpha, beta, a, b);
  }
  return lp;
}

//with weight
// [[Rcpp::export]]
List doVB_pois_s_2D_ww(const arma::vec & y,
                    const arma::uvec & rowi,
                    const arma::uvec & coli,
                    const int & L,
                    const int & iter,
                    const int & subiter,
                    const arma::field<arma::vec> & weight,
                    const double & a,
                    const double & b,
                    const double & N1,
                    const int & Nr,
                    const int & Nc,
                    const int & bsize,
                    const arma::vec & lr_param,
                    const std::string & lr_type,
                    const bool & display_progress){
  arma::field<arma::mat> logV(2);
  arma::field<arma::mat> alpha(2);
  arma::field<arma::mat> V(2);
  V(0) = arma::randg<arma::mat>(Nr,L);
  V(1) = arma::randg<arma::mat>(Nc,L);
  logV(0) = log(V(0));
  logV(1) = log(V(1));
  alpha(0) = arma::ones<arma::mat>(Nr, L);
  alpha(1) = arma::ones<arma::mat>(Nc, L);
  arma::mat beta(2,L);
  beta.row(0) = sum(V(1), 0);
  beta.row(1) = sum(V(0), 0);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  const double NS = N1 / ((double) bsize);
  double invS = 1.0 / ( (double) bsize ); 
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
    //rho *= invS;
    const double NS = N1 / ((double) bsize);
    double invS = 1.0 / ( (double) bsize ); 
    for(int step = 0; step < bags.n_cols; step++){
      arma::uvec bag = sort(bags.col(step));
      arma::vec Sy = y.rows(bag);
      arma::umat SX(bag.n_rows, 2);
      SX.col(0) = rowi.rows(bag);
      SX.col(1) = coli.rows(bag);
      arma::field<arma::uvec> uid(2);
      uid(0) = unique(SX.col(0));
      uid(1) = unique(SX.col(1));
      arma::field<arma::mat> alpha_s = alpha;
      arma::mat beta_s = beta;
      lp(epoc) += doVB_pois_s_sub_2D(Sy, SX, L, subiter, a, b, NS, uid, weight, alpha_s, beta_s, V, logV);
      up_vpar(rho, rho2, uid, alpha, alpha_s, beta, beta_s);
    }
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}
