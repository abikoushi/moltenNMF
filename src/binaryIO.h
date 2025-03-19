#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

void readRowsFromBinary_umat(arma::umat & X,
                             const std::string & filepath,
                             const arma::uvec & rows);

void readRowsFromBinary(arma::mat & Y, 
                        const std::string & filepath, 
                        const arma::uvec & rows);
