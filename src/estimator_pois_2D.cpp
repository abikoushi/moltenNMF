#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "logexpfuns.h"
#include "rand.h"
#include "up_shape_2D.h"
#include "lr.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

double up_B_2D(const arma::field<arma::mat> & alpha,
               arma::mat & beta,
               arma::field<arma::mat> & V,
               arma::field<arma::mat> & logV,
               const double & b,
               const int & L,
               const arma::uvec & dims, 
               const int & k){
  double lp = 0;
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  for(int l=0; l<L; l++){
    double B0 = sum(V(not_k).col(l));
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
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 0);
  lp += up_B_2D(alpha, beta, V, logV, b, L, dims, 0);
  lp += up_A_2D(alpha, logV, y, X, a, L, dims, 1);
  lp += up_B_2D(alpha, beta, V, logV, b, L,dims, 1);
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
    double lp0 = up_theta_2D(alpha, beta, V, logV, y, X, L, dims, a, b);
    lp(i) = lp0 + kld2(alpha, beta, a, b);
    pb.increment();
  }
  lp -= sum(lgamma(y+1));
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}

////
//SVB
///

double up_Bs_2D(const arma::field<arma::mat> & alpha,
                arma::mat & beta,
                arma::field<arma::mat> & V,
                arma::field<arma::mat> & logV,
                arma::mat & SumV,
                const arma::field<arma::uvec> & uid,
                const double & b,
                const int & L,
                const double & NS,
                const int & k){
  double lp = 0;
  int not_k = 1;
  if(k==1){
    not_k = 0;
  }
  for(int l=0; l<L; l++){
    arma::vec vl = V(not_k).col(l);
    vl = vl.rows(uid(not_k));
    SumV(not_k,l) += NS*sum(vl);
    //double B0 = NS*(sum(V(not_k).col(l)));
    lp -= SumV(not_k,l) - b;
    beta(k,l) = SumV(not_k,l);
    //beta(k,l) = B0 + b;
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


double up_theta_s_2D(arma::field<arma::mat> & alpha,
                     arma::mat & beta,
                     arma::field<arma::mat> & V,
                     arma::field<arma::mat> & logV,
                     arma::mat & SumV,
                     const int & L,
                     const arma::vec & y,
                     const arma::umat & X,
                     const arma::field<arma::uvec> & uid,
                     const double & a,
                     const double & b,
                     const double & NS){
  double lp = 0;
  //rows
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 0);
  lp += up_Bs_2D(alpha, beta, V, logV, SumV, uid, b, L, NS, 0);
  //columns
  lp += up_As_2D(alpha, logV, y, X, a, L, uid, 1);
  lp += up_Bs_2D(alpha, beta, V, logV, SumV, uid, b, L, NS,  1);
  return lp;
}

//use SumV
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
  arma::mat SumV(2,L);
  for(int j=0;j<2;j++){
    SumV.row(j) = (beta.row(j)) - sum(V(j).rows(uid(j)), 0)*NS;
  }
  SumV.elem(find(SumV<0.0)).fill(0.0);
  double lp = 0.0;
  for (int i=0; i<iter; i++) {
    double lp1 = up_theta_s_2D(alpha, beta, V, logV, SumV, L, y, X, uid, a, b, NS);
    lp = lp1 + kld2(alpha, beta, a, b);
  }
  return lp;
}


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
  arma::mat beta = arma::ones<arma::mat>(2,L);
  arma::vec lp = arma::zeros<arma::vec>(iter);
  std::unique_ptr<lr> g;
  if(lr_type == "const"){
    g.reset(new lr_const);
  }else if(lr_type == "power"){
    g.reset(new lr_power);
  }else{
    Rcpp::stop("This lr_type is not implemented\n");
  }
  const double NS = (bsize/N1);
  Progress pb(iter, display_progress);
  for(int epoc=0; epoc<iter; epoc++){
    arma::umat bags = randpick_c(N1, bsize);
    double rho = g -> lr_t(epoc, lr_param);
    double rho2 = 1.0 - rho;
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
      lp(epoc) += doVB_pois_s_sub_2D(Sy, SX, L, subiter, a, b, NS, uid, alpha_s, beta_s, V, logV);
      for(int j = 0; j < 2; j++){
        alpha(j).rows(uid(j)) = alpha(j).rows(uid(j))*rho2 + rho*alpha_s(j).rows(uid(j));
      }
      beta = rho2*beta + rho*beta_s;
    }
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=lp);
}
