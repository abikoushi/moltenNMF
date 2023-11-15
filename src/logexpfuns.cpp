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