#include <RcppArmadillo.h>
#include "subset_sp.h"

void subset_spx(arma::uvec & xi2, arma::uvec & xp2,
                const arma::uvec & xi, const arma::uvec & xp,
                const arma::uvec & uid) {
  xp2 = xp;
  std::vector<arma::uword> result;
  int count_p = 0;
  for(int i = 0; i < ((int) xp.n_rows) -1; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      if(any(xi(j) == uid)){
        count_p += 1;
        arma::uvec newxi = find(xi(j) == uid);
        result.push_back(newxi(0));
      }
    }
    xp2(i+1) = count_p;
  }
  xi2 = arma::conv_to<arma::uvec>::from(result);
}

arma::uvec subset_xp(const arma::uvec & xi, const arma::uvec & xp,
                     const arma::uvec & uid) {
  arma::uvec xp2 = xp;
  int count_p = 0;
  for(int i = 0; i < ((int) xp.n_rows) - 1; i++){
    for(int j = xp(i); j < (int) xp(i+1); j++){
      if(any(xi(j) == uid)){
        count_p += 1;
      }
    }
    xp2(i+1) = count_p;
  }
  return xp2;
}


void filter_y(arma::vec & yv2, arma::uvec & yi2,
              const arma::vec & yv, const arma::uvec & yi,
              const arma::uvec & ind) {
  std::vector<double> result_v;
  std::vector<int> result_i;
  int count_i = 0; 
  for(int i = 0; i < (int) ind.n_rows; i++){
    if(any(yi==ind(i))){
      arma::uvec ti = find(yi==ind(i));
      result_v.push_back(yv(ti(0)));
      result_i.push_back(count_i);
    }
    count_i++;
  }
  yv2 = arma::conv_to<arma::vec>::from(result_v);
  yi2 = arma::conv_to<arma::uvec>::from(result_i);
  //return Rcpp::List::create(yi2,yv2);
}
