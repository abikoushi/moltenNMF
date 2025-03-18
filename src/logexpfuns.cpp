#include <RcppArmadillo.h>
#include "logexpfuns.h"

double logsumexp(const arma::rowvec & x){
  double maxx = max(x);
  return maxx + std::log(sum(exp(x-maxx)));
}

arma::rowvec softmax(const arma::rowvec & x){
  double maxx = max(x);
  arma::rowvec u = x-maxx;
  return exp(u)/sum(exp(u));
}

arma::mat mat_digamma(arma::mat a){
  int K = a.n_rows;
  int L = a.n_cols;
  arma::mat out(K,L);
  for(int k=0;k<K;k++){
    for(int l=0;l<L;l++){
      out(k,l) = R::digamma(a(k,l));
    }
  }
  return out;
}

arma::vec vec_digamma(arma::vec a){
  int K = a.n_rows;
  arma::vec out(K);
  for(int k=0;k<K;k++){
    out[k] = R::digamma(a[k]);
  }
  return out;
}

void up_log_gamma(arma::mat & logv, const arma::vec & a, const double & logb, const int & l){
  int K = logv.n_rows;
  for(int k=0; k<K; k++){
    logv(k,l) = R::digamma(a(k)) - logb;
  }
}

double klgamma_sub(double a, double b, double c, double d){
  return - c * d / a - b * log(a) - lgamma(1/b) + (b-1)*(R::digamma(1/d) + log(c));
}

double kl2gamma(double a1, double b1, double a2, double b2){
  return klgamma_sub(a2,b2,a2,b2) - klgamma_sub(a1,b1,a2,b2);
}

double kld(const arma::mat & alpha,
           const arma::rowvec & beta,
           const double & a,
           const double & b){
  double lp = 0;
  for(int i=0; i<alpha.n_rows; i++){
    for(int l=0; l<alpha.n_cols; l++){
      lp += kl2gamma(a, b, alpha(i,l), beta(l));
    }
  }
  return lp;
}

double kld2(const arma::field<arma::mat> & alpha,
            const arma::mat & beta,
            const double & a,
            const double & b){
  double lp = 0;
  int K=alpha.n_rows;
  for(int k=0; k<K ;k++){
    arma::mat alpha_k = alpha(k);
    arma::rowvec beta_k = beta.row(k);
    lp += kld(alpha_k, beta_k, a, b);    
  }
  return lp;
}
