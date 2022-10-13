#include <RcppArmadillo.h>
#include "myproduct.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

double logsumexp(const arma::rowvec & x){
  double maxx = max(x);
  return maxx + std::log(sum(exp(x-maxx)));
}

void up_x(arma::mat & Xprob,
           const arma::vec & y,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & miss_row,
           const arma::uvec & miss_col,
           const arma::uvec & varind,
           const arma::uvec & dropK,
           const arma::mat & loglambda,
           const arma::mat & lambda,
           const int & D,
           const int & L,
           const arma::vec & etahat){
  arma::mat et0 = myprod(y.n_rows, xi, xp, lambda);
  arma::mat t0 = mysum(y.n_rows, xi, xp, loglambda);
  double den = R::digamma(sum(etahat));
  arma::vec logw = vec_digamma(etahat) - den;
  //arma::vec log1mw = vec_digamma(etahat) - den;
  for(int i=0; i<miss_row.n_rows; i++){
    arma::mat vt = lambda.rows(varind[miss_col[i]], varind[miss_col[i]+1]-1);
    arma::mat lvt = loglambda.rows(varind[miss_col[i]], varind[miss_col[i]+1]-1);
    double den = dropK[varind[miss_col[i]]];
    for(int j=varind[miss_col[i]]; j<varind[miss_col[i]+1]-1; j++){
      arma::rowvec et1 = vt.row(j-varind[miss_col[i]]) % et0.row(miss_row[i]);
      arma::rowvec t1 = lvt.row(j-varind[miss_col[i]]) + t0.row(miss_row[i]);
      double T1 = arma::as_scalar(y.row(miss_row[i]))*logsumexp(t1)-arma::as_scalar(sum(et1))+logw[j];
      //double T0 = arma::as_scalar(y.row(miss_row[i]))*logsumexp(t0.row(miss_row[i]))-arma::as_scalar(sum(et0.row(miss_row[i])))+log1mw[j];
      den += exp(T1);
      Xprob.row(i).col(j) = exp(T1);
    }
    Xprob.row(i) /= den;
  }
}

// [[Rcpp::export]]
List doVB_pois_missing(const arma::vec & y,
                       const arma::uvec & xi,
                       const arma::uvec & xp,
                       const arma::uvec & miss_row,
                       const arma::uvec & miss_col,
                       const arma::vec & sumx,
                       const arma::uvec & varind,
                       const arma::uvec & dropK,
                       const int & D,
                       const int & L,
                       const int & iter,
                       const double & a,
                       const double & b,
                       const double & eta){
  int N = y.n_rows;
  int K = varind.n_rows - 1;
  arma::mat lambda = arma::randg<arma::mat>(D,L);
  arma::vec etahat = sumx+eta;
  //arma::vec etahat2 = sum1mx+eta;
  arma::mat Xprob = arma::zeros<arma::mat>(miss_row.n_rows, D);
  arma::mat loglambda = log(lambda);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec ll(iter);
  for (int i=0; i<iter; i++) {
    up_x(Xprob, y, xi, xp, miss_row, miss_col, varind, dropK, loglambda, lambda, D, L, etahat);
    etahat = sumx + sum(Xprob,0).t();
    //etahat2 = sum1mx + sum(1.0-Xprob,0).t();
    arma::mat r =  myprod(N, xi, xp, exp(loglambda));
    for (int j=0; j<miss_row.n_rows; j++){
      r.row(miss_row[j]) %=  exp(Xprob.row(j)*loglambda); 
    }
    arma::vec den = sum(r, 1);
    alpha = mysum_t(D, xi, xp, r.each_col()%(y/den)) + a;
    for(int k=0; k<K; k++){
      arma::mat tmpXl = myprod_skip(N, xi, xp, lambda, varind[k], varind[k+1]);
      arma::mat B = mysum_t(varind[k+1] - varind[k], xi, xp.rows(varind[k], varind[k+1]), tmpXl) + b;
      beta.rows(varind[k], varind[k+1]-1) = B;
      lambda.rows(varind[k], varind[k+1]-1) = alpha.rows(varind[k],varind[k+1]-1)/B;
    }
    loglambda = mat_digamma(alpha) - log(beta);
    ll.row(i) = sum(-den + y%log(den) - lgamma(y+1)) +
      + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
      - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha));
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("Xprob")=Xprob,
                      Named("shape_x")=etahat,
                      Named("ELBO")=ll);
}