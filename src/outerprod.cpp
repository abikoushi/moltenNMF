#include "RcppArmadillo.h"
#include "outerprod.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double sumouter2vec(arma::vec & x, arma::vec & y){
  double res = 0;
  for(int i = 0; i < (int) x.n_rows; i++){
    res += sum(x(i) * y);
  }
  return res;
}

double sumouterprod_s(const arma::field<arma::mat> & V,
                      const int & k,
                      const int & l,
                      const arma::field<arma::uvec> & uid){
  int K = V.n_rows;
  uvec dummy = linspace<uvec>(0, K-1, K);
  dummy.shed_row(k);
  arma::vec V_jl = V(dummy(0)).col(l);
  arma::vec tout = V_jl.rows(uid(dummy(0)));
  for(int j=1; j < ((int) dummy.n_rows) - 1; j++){
    arma::uvec uid_j = uid(dummy(j));
    V_jl = V(dummy(j)).col(l);
    arma::vec V1 = V_jl.rows(uid_j);
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows > 1){
    V_jl = V(dummy(dummy.n_rows-1)).col(l);
    arma::vec VK = V_jl.rows(uid(dummy(dummy.n_rows-1)));
    out = sumouter2vec(tout, VK);
  }else{
    out = sum(tout);
  }
  return out;
}

double sumouterprod_s(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l,
                    const arma::field<arma::vec> & weight,
                    const arma::field<arma::uvec> & uid){
  int K = V.n_rows;
  uvec dummy = linspace<uvec>(0, K-1, K);
  dummy.shed_row(k);
  arma::vec weight0 = weight(dummy(0));
  arma::vec V_jl = V(dummy(0)).col(l);
  arma::vec f_j =  weight(dummy(0)).rows(uid(dummy(0)));
  arma::vec tout = V_jl.rows(uid(dummy(0)))%f_j;
  for(int j = 1; j < ((int)dummy.n_rows) - 1; j++){
    f_j = weight(dummy(j)).rows(uid(dummy(j)));
    V_jl = V(dummy(j)).col(l);
    arma::vec V1 = V_jl.rows(uid(dummy(j))) % f_j;
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows > 1){
    arma::vec weightK = weight(dummy(dummy.n_rows-1)).rows(uid(dummy(dummy.n_rows-1)));
    V_jl = V(dummy(dummy.n_rows-1)).col(l);
    arma::vec VK = V_jl.rows(uid(dummy(dummy.n_rows-1))) % weightK;
    out = sumouter2vec(tout, VK);
  }else{
    out = sum(tout);
  }
  return out;
}

double sumouterprod(const arma::field<arma::mat> & V,
                    const int & k,
                    const int & l){
  int K = V.n_rows;
  uvec dummy = linspace<uvec>(0, K-1, K);
  dummy.shed_row(k);
  arma::vec tout = V(dummy(0)).col(l);
  for(int j = 1; j < ((int) dummy.n_rows) - 1; j++){
    arma::vec V1 = V(dummy(j)).col(l);
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows > 1){
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
  for(int j=1; j < ((int) dummy.n_rows) - 1; j++){
    arma::vec weight1 = weight(dummy(0));
    arma::vec V1 = V(dummy(j)).col(l) % weight1;
    tout = vectorise(tout * V1.t());
  }
  double out;
  if(dummy.n_rows > 1){
    arma::vec weightK = weight(dummy(dummy.n_rows-1));
    arma::vec VK = V(dummy(dummy.n_rows-1)).col(l) % weightK;
    out = sumouter2vec(tout, VK);
  }else{
    out = sum(tout);
  }
  return out;
}

double sumouterprod_mat(const arma::mat & V,
                        const arma::uvec & varind,
                        const int & k,
                        const int & l,
                        const arma::field<arma::vec> & weight){
  int K = varind.n_rows-1;
  arma::vec tout = arma::ones<arma::vec>(1);
  for(int j = 0; j < k; j++){
    arma::vec weight1 = weight(j);
    arma::vec Vkl = V.rows(varind(j), varind(j+1)-1).col(l) % weight1;
    tout = vectorise(tout * Vkl.t());
    //Rprintf("%d ", j);
  }
  //skip j==k
  for(int j = k + 1; j < K; j++){
    arma::vec weight1 = weight(j);
    arma::vec Vkl = V.rows(varind(j), varind(j+1)-1).col(l) % weight1;
    tout = vectorise(tout * Vkl.t());
    //Rprintf("%d ", j);
  }
  //Rprintf("\n");
  return sum(tout);
}

double sumouterprod_mat(const arma::mat & V,
                        const arma::uvec varind,
                        const int & k,
                        const int & l,
                        const arma::field<arma::vec> & prob,
                        const arma::vec & freq){
  int K = varind.n_rows-1;
  arma::vec tout = arma::ones<arma::vec>(1);
  // double logN_k = 0.0;
  for(int j = 0; j < k; j++){
    arma::vec weight1 = exp(log(freq(j))+log(prob(j)));
    arma::vec Vkl = V.rows(varind(j), varind(j+1)-1).col(l) % weight1;
    tout = vectorise(tout * Vkl.t());
  }
  for(int j = k + 1; j < K; j++){
    arma::vec weight1 = exp(log(freq(j))+log(prob(j)));
    arma::vec Vkl = V.rows(varind(j), varind(j+1)-1).col(l) % weight1;
    tout = vectorise(tout * Vkl.t());
  }
  return sum(tout);
}
