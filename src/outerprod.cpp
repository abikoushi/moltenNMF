#include "RcppArmadillo.h"
#include "outerprod.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double sumouter2vec(arma::vec & x, arma::vec & y){
  double res = 0;
  for(int i = 0; i < x.n_rows; i++) {
    res += sum(x(i) * y);
  }
  return res;
}

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l){
  int K = V.n_rows;
  uvec dummy = linspace<uvec>(0, K-1, K);
  dummy.shed_row(k);
  arma::vec tout = V(dummy(0)).col(l);
  for(int j=1; j<dummy.n_rows-1; j++){
    arma::vec V1 = V(dummy(j)).col(l);
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows>1){
    arma::vec VK = V(dummy(dummy.n_rows-1)).col(l);
    out = sumouter2vec(tout, VK);
  }else{
    out = sum(tout);
  }
  return out;
}

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l,
                    const arma::field<arma::vec> & weight){
  int K = V.n_rows;
  uvec dummy = linspace<uvec>(0, K-1, K);
  dummy.shed_row(k);
  arma::vec weight0 = weight(dummy(0));
  arma::vec tout = V(dummy(0)).col(l) % weight0;
  for(int j=1; j<dummy.n_rows-1; j++){
    arma::vec weight1 = weight(dummy(0));
    arma::vec V1 = V(dummy(j)).col(l) % weight1;
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows>1){
    arma::vec weightK = weight(dummy(dummy.n_rows-1));
    arma::vec VK = V(dummy(dummy.n_rows-1)).col(l) % weightK;
    out = sumouter2vec(tout, VK);
  }else{
    out = sum(tout);
  }
  return out;
}