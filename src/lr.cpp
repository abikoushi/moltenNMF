#include "lr.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

void set_lr_method(std::unique_ptr<lr> & g, const std::string & lr_type){
  if(lr_type == "const"){
    g.reset(new lr_const);
  }else if(lr_type == "exponential"){
    g.reset(new lr_power);
  }else{
    Rcpp::stop("This lr_type is not implemented\n");
  }
}

// [[Rcpp::export]]
double check_lr(const int & epoc,
                const arma::vec & lr_param,
                const std::string & lr_type){
  std::unique_ptr<lr> g;
  set_lr_method(g, lr_type);
  double rho = g -> lr_t(epoc, lr_param);
  return rho;
}
