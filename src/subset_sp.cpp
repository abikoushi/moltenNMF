#include <RcppArmadillo.h>
#include "subset_sp.h"

void subset_spx(arma::uvec & xi2, arma::uvec & xp2,
                const arma::uvec & xi, const arma::uvec & xp,
                const arma::uvec & uid, arma::uvec & up) {
  std::unordered_map<arma::uword, arma::uword> uid_lookup;
  for (arma::uword i = 0; i < uid.n_elem; ++i) {
    uid_lookup[uid(i)] = i;
  }
  
  std::vector<arma::uword> xi_buf;
  std::vector<arma::uword> up_buf;
  
  xp2 = arma::uvec(xp.n_rows, arma::fill::zeros);
  arma::uword count = 0;
  
  for (arma::uword i = 0; i < xp.n_rows - 1; ++i) {
    for (arma::uword j = xp(i); j < xp(i + 1); ++j) {
      arma::uword val = xi(j);
      auto it = uid_lookup.find(val);
      if (it != uid_lookup.end()) {
        xi_buf.push_back(it->second);
        up_buf.push_back(i);
        ++count;
      }
    }
    xp2(i + 1) = count;
  }
  
  xi2 = arma::uvec(xi_buf);
  up = arma::unique(arma::uvec(up_buf));
}


// void subset_spx_old(arma::uvec & xi2, arma::uvec & xp2,
//                 const arma::uvec & xi, const arma::uvec & xp,
//                 const arma::uvec & uid, arma::uvec & up) {
//   xp2 = xp;
//   std::vector<arma::uword> vector_up;
//   std::vector<arma::uword> result;
//   int count_p = 0;
//   for(arma::uword i = 0; i < (xp.n_rows - 1); i++){
//     for(arma::uword j = xp(i); j < xp(i+1); j++){
//       if(any(xi(j) == uid)){
//         count_p += 1;
//         arma::uvec newxi = find(xi(j) == uid);
//         result.push_back(newxi(0));
//         vector_up.push_back(i);
//       }
//     }
//     xp2(i+1) = count_p;
//   }
//   xi2 = arma::conv_to<arma::uvec>::from(result);
//   up = unique(arma::conv_to<arma::uvec>::from(vector_up));
// }

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
