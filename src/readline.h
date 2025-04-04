#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void readmtx(arma::umat & X,
             arma::vec & val,
             const std::string & readtxt,
             const arma::uvec & bag,
             const int & n_header);

void rankindex(arma::uvec & x, const arma::uvec & uid);
