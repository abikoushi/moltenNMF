#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

class lr {
public:
  virtual double lr_t(const int t, const arma::vec & lrparam)=0;
  virtual ~lr(){}
};

class lr_power : public lr{
  double lr_t(const int t, const arma::vec & lrparam){
    return pow(t+lrparam(0), -lrparam(1));
  }
};

class lr_const : public lr{
  double lr_t(const int t, const arma::vec & lrparam){
    return lrparam(0);
  }
};