#include "RcppArmadillo.h"
#include "up_shape_2D.h"
// [[Rcpp::depends(RcppArmadillo)]]

double up_A_2D(arma::field<arma::mat> & alpha,
            const arma::field<arma::mat> & logV,
            const arma::vec & y,
            const arma::umat & X,
            const double & a,
            const int & L,
            const arma::uvec & dims,
            const int & k){
  //initialize by hyper parameter
  arma::mat alpha_k(dims(k),L);
  alpha_k.fill(a);
  alpha(k) = alpha_k;
  //inclement sufficient statistics
    arma::mat r = arma::zeros<arma::mat>(y.n_rows, L);
    for(int j = 0; j < (int) X.n_cols; j++){
      arma::mat logVj = logV(j);
      r += logVj.rows(X.col(j));
    }
    r = exp(r);
    arma::vec R = sum(r,1);
    r.each_col() /= R;
    r.each_col() %= y;
  alpha(k).rows(X.col(k)) +=  r;
  return sum(y%log(R)-R);
}

double up_A_2D(arma::field<arma::mat> & alpha,
               const arma::field<arma::mat> & logV,
               const arma::vec & y,
               const arma::umat & X,
               const double & a,
               const int & L,
               const arma::uvec & dims){
  int K = X.n_cols;
  //calculate latent statistics
  arma::mat r = arma::zeros<arma::mat>(y.n_rows, L);
  for(int j = 0; j < K; j++){
    arma::mat logVj = logV(j);
    r += logVj.rows(X.col(j));
  }
  r = exp(r);
  arma::vec R = sum(r, 1);
  r.each_col() /= R;
  r.each_col() %= y;
  for(int k = 0; k < K; k++){
    //initialize by hyper parameter
    arma::mat alpha_k(dims(k),L);
    alpha_k.fill(a);
    alpha(k) = alpha_k;
    //inclement sufficient statistics
    alpha(k).rows(X.col(k)) +=  r;
  }
  return sum(y%log(R)-R);
}

////
//Use stochastic mini-batches
////

double up_As_2D(arma::field<arma::mat> & alpha,
                const arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const double & a,
                const int & L,
                const arma::field<arma::uvec> & uid,
                const int & k,
                const double & NS){
  //initialize by hyper parameter
  arma::uvec uid_k = uid(k);
  arma::mat alpha_k = alpha(k);
  alpha_k.rows(uid_k).fill(a);
  alpha(k) = alpha_k;
  //inclement sufficient statistics
  arma::mat r = arma::zeros<arma::mat>(y.n_rows, L);
  for(int j = 0; j < (int) X.n_cols; j++){
    arma::mat logVj = logV(j);
    r += logVj.rows(X.col(j));
  }
  r = exp(r);
  arma::vec R = sum(r, 1);
  r.each_col() /= R;
  r.each_col() %= y;
  alpha(k).rows(X.col(k)) += r * NS;
  return sum(y % log(R)-R); //lp
}

double up_As_2D(arma::field<arma::mat> & alpha,
                const arma::field<arma::mat> & logV,
                const arma::vec & y,
                const arma::umat & X,
                const double & a,
                const int & L,
                const arma::field<arma::uvec> & uid,
                const double & NS){
  int K= X.n_cols;
  //calculate latent statistics
  arma::mat r = arma::zeros<arma::mat>(y.n_rows, L);
    for(int j = 0; j < K; j++){
      arma::mat logVj = logV(j);
      r += logVj.rows(X.col(j));
    }
    r = exp(r);
    arma::vec R = sum(r, 1);
    r.each_col() /= R;
    r.each_col() %= y;
  for(int k = 0; k < K; k++){
      arma::uvec uid_k = uid(k);
      arma::mat alpha_k = alpha(k);
      //initialize by hyper parameter
      alpha_k.rows(uid_k).fill(a);
      alpha(k) = alpha_k;
      //inclement sufficient statistics
      alpha(k).rows(X.col(k)) += r * NS;    
  }
  return  sum(y % log(R)); //lp
}

//////
//The same applies for multidimensional arrays.
//////
