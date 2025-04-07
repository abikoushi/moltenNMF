#include "readline.h"
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <string>
#include <vector>
using namespace Rcpp;

void readmtx(arma::umat & X,
             arma::vec & val,
             const std::string & readtxt,
             const arma::uvec & bag,
             const int & n_header) {
  std::ifstream file(readtxt);
  std::string str;    
  int index = 0;
  int n = 0;
  while (std::getline(file, str)){
      if(index == ((int) bag(n)) + n_header){
        std::stringstream ss(str);
        std::vector<std::string> svec;
        while( ss.good() ){
          std::string substr;
          getline(ss, substr, ' ');
          svec.push_back(substr);
        }
        int dim = svec.size();
        for(int k = 0; k < dim-1; k++){
          int x;
          x = stoi(svec[k]);
          X(n,k) = x-1;
        }
        double v;
        v = stod(svec[dim-1]);
        val(n) = v;
        n++;
      }
      index++;
      if(n >= (int) bag.n_rows){
        break;
      }
  }
}

// [[Rcpp::export]]
List read_mtx(const std::string & readtxt,
              const arma::uvec & bag,
              const int & x_dim){
  arma::umat X(bag.n_rows, x_dim);
  arma::vec val(bag.n_rows);
  readmtx(X, val, readtxt, bag, 2);
  return List::create(X, val);
}

void rankindex(arma::uvec & x, const arma::uvec & uid){
  for(int i=0; i < (int) uid.n_rows; i++){
    arma::uvec idx = find(uid(i) == x);
    x.rows(idx).fill(i);
  }
}
