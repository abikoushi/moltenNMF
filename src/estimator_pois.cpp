#include <RcppArmadillo.h>
#include "myproduct.h"
#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double logsumexp(const arma::rowvec & x){
  double maxx = max(x);
  return maxx + std::log(sum(exp(x-maxx)));
}

arma::rowvec softmax(const arma::rowvec & x){
  double maxx = max(x);
  arma::rowvec u = x-maxx;
  return exp(u)/sum(exp(u));
}

///
void up_A_r(arma::mat & alpha,
            arma::vec & R,
            const arma::mat & loglambda,
            const arma::vec & y,
            const arma::uvec & xj,
            const arma::uvec & xp,
            const double & a){
  arma::mat r =  myprod_r(y.n_rows, xj, xp, exp(loglambda)); //N,L
  R = sum(r, 1);
  alpha = mysum_t_r(alpha.n_rows, xj, xp, r.each_col()%(y/R)) + a; //D,L
}

void up_B_r(const int & N,
            arma::mat & beta,
            arma::mat & V,
            arma::mat & logV,
            const arma::mat & alpha,
            const arma::uvec & xj,
            const arma::uvec & xp,
            const arma::uvec & varind,
            const double & b){
  int K = varind.n_rows - 1;
  for(int k=0; k<K; k++){
    arma::mat tmpXl = myprod_skip(N, xj, xp, V, varind[k], varind[k+1]);
    arma::mat B = mysum_t_r(varind[k+1]-varind[k], xj, xp.rows(varind[k], varind[k+1]), tmpXl) + b;
    beta.rows(varind[k], varind[k+1]-1) = B;
    V.rows(varind[k], varind[k+1]-1) = alpha.rows(varind[k],varind[k+1]-1)/B;
  }
  logV = mat_digamma(alpha) - log(beta); 
}
///

void up_B2(const int & N,
           arma::mat & beta,
           arma::mat & V,
           arma::mat & logV,
           arma::vec & z,
           const arma::mat & alpha,
           const arma::uvec & xi,
           const arma::uvec & xp,
           const arma::uvec & varind,
           const double & b){
  int K = varind.n_rows - 1;
  for(int k=0; k<K; k++){
    arma::mat tmpXl = myprod_skip(N, xi, xp, V, varind[k], varind[k+1]);
    tmpXl.each_col() %= z;
    arma::mat B = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k], varind[k+1]), tmpXl) + b;
    beta.rows(varind[k], varind[k+1]-1) = B;
    V.rows(varind[k], varind[k+1]-1) = alpha.rows(varind[k],varind[k+1]-1)/B;
  }
  logV = mat_digamma(alpha) - log(beta); 
}

//shape parameters
void up_A(arma::mat & alpha,
          arma::vec & R,
          const arma::mat & loglambda,
          const arma::vec & y,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const double & a){
  arma::mat r =  myprod(y.n_rows, xi, xp, exp(loglambda)); //N,L
  R = sum(r, 1);
  alpha = mysum_t(alpha.n_rows, xi, xp, r.each_col()%(y/R)) + a; //D,L
}

//rate parameters
void up_B(const int & N,
          arma::mat & beta,
          arma::mat & V,
          arma::mat & logV,
          const arma::mat & alpha,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::uvec & varind,
          const double & b){
  int K = varind.n_rows - 1;
  int L = V.n_cols;
  for(int l=0;l<L;l++){
    arma::vec vl = myprodvec(N, xi, xp, V.col(l));
    for(int k=0; k<K; k++){
      vl /= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
      arma::vec B = mysum_t(varind[k+1]-varind[k], xi, xp.rows(varind[k], varind[k+1]), vl) + b;
      beta.col(l).rows(varind[k], varind[k+1]-1) = B;
      V.col(l).rows(varind[k], varind[k+1]-1) = alpha.col(l).rows(varind[k],varind[k+1]-1)/B;
      vl %= myprodvec_sub(N, xi, xp, varind[k], varind[k+1], V.col(l));
    }
    logV = mat_digamma(alpha) - log(beta);
  }
}

//precision parameter
void up_tau(double & tau,
            arma::vec & z,
            const arma::vec & y,
            const arma::mat & V,
            const arma::uvec & xi,
            const arma::uvec & xp){
  arma::vec ybar = sum(myprod(z.n_rows, xi, xp, V), 1);
  arma::vec ahat = (y+tau);
  arma::vec bhat = (ybar+tau);
  arma::vec logz = vec_digamma(ahat) - log(bhat);
  z = ahat/bhat;
  tau += sum(logz - z + log(tau) + 1 - R::digamma(tau))/(y.n_rows*(1/tau - R::trigamma(tau)));
}

void up_x(arma::mat & Xprob,
          const arma::vec & y,
          const arma::uvec & xi,
          const arma::uvec & xp,
          const arma::vec & sumx,
          const arma::uvec & miss_row,
          const arma::uvec & miss_col,
          const arma::uvec & varind,
          const arma::mat & loglambda,
          const arma::mat & lambda,
          const int & D,
          const int & L,
          arma::vec & logw,
          arma::vec & etahat,
          const double & eta){
  Xprob.fill(0);
  arma::mat et0 = myprod(y.n_rows, xi, xp, lambda);
  arma::mat t0 = mysum(y.n_rows, xi, xp, loglambda);
  double den = R::digamma(sum(etahat));
  logw = vec_digamma(etahat) - den;
  for(int i=0; i<miss_row.n_rows; i++){
    //Rprintf("%d: ", i);
    arma::mat vt = lambda.rows(varind[miss_col[i]], varind[miss_col[i]+1]-1);
    arma::mat lvt = loglambda.rows(varind[miss_col[i]], varind[miss_col[i]+1]-1);
    //double den = 0;
    int Jmax = varind[miss_col[i]+1]-varind[miss_col[i]]-1;
    arma::rowvec lp = arma::zeros<arma::rowvec>(Jmax+1); 
    for(int j=0; j<Jmax; j++){
      arma::rowvec et1 = vt.row(j) % et0.row(miss_row[i]);
      arma::rowvec t1 = lvt.row(j) + t0.row(miss_row[i]);
      lp(j) = arma::as_scalar(y.row(miss_row[i]))*logsumexp(t1)-arma::as_scalar(sum(et1))+logw[varind[miss_col[i]]+j];
    }
    Xprob.row(i).cols(varind[miss_col[i]], varind[miss_col[i]+1]-1) = softmax(lp);
  }
  etahat = sumx + sum(Xprob, 0).t() + eta;
}

double lowerbound_logML_pois(const arma::mat & alpha,
                        const arma::mat & beta,
                        const arma::mat & lambda,
                        const arma::mat & loglambda,
                        const arma::vec & R,
                        const arma::vec & y,
                        const double & a,
                        const double & b){
  return sum(-R + y%log(R) - lgamma(y+1)) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha));
}

double lowerbound_logML_mult(const arma::vec & logW,
                             const arma::vec & sumx,
                             const arma::mat & Xprob,
                             const arma::vec & etahat,
                             const double & eta){
  return sum((sumx+arma::trans(sum(Xprob,0)))%logW) + sum(logW*eta) -
    (lgamma(logW.n_cols*eta) - logW.n_cols*lgamma(eta) +
    - lgamma(sum(etahat)) + sum(lgamma(etahat)) +
    sum((eta - etahat)%logW));
}

double lowerbound_logML2(const arma::vec & z,
                         const arma::mat & alpha,
                        const arma::mat & beta,
                        const double & tau,
                        const arma::mat & lambda,
                        const arma::mat & loglambda,
                        const arma::vec & R,
                        const arma::vec & y,
                        const arma::uvec & xi,
                        const arma::uvec & xp,
                        const double & a,
                        const double & b){
  arma::vec ahat = (y+tau);
  arma::vec ybar = sum(myprod(y.n_rows, xi, xp, lambda), 1);
  arma::vec bhat = (ybar+tau);
  arma::vec logz = vec_digamma(ahat) - log(bhat);
  return sum(-R + y%log(R) - lgamma(y+1) + y%logz) +
    + accu((a-1)*loglambda - b*lambda + a*log(beta) - std::lgamma(a)) +
    - accu((alpha-1)%loglambda - beta%lambda + alpha%log(beta) - lgamma(alpha))+
    + sum(z%bhat - z*a - ahat%log(bhat) + a*log(a)+lgamma(ahat)-lgamma(a));;
}

// [[Rcpp::export]]
List doVB_pois(const arma::vec & y,
               const arma::uvec & xi,
               const arma::uvec & xp,
               const arma::uvec & varind,
               const int & D,
               const int & L,
               const int & iter,
               const double & a,
               const double & b,
               const bool & display_progress){
  int N = y.n_rows;
  arma::mat V = arma::randg<arma::mat>(D,L);
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  arma::vec ll(iter);
  Progress pb(iter, display_progress);
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    up_B(N, beta, V, logV, alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, V, logV, R, y, a, b);
    pb.increment();
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("ELBO")=ll);
}

// [[Rcpp::export]]
List doVB_negbin(arma::vec y,
                 arma::uvec xi,
                 arma::uvec xp,
                 arma::uvec varind,
                 int D,
                 int L,
                 int iter,
                 double a,
                 double b){
  int N = y.n_rows;
  arma::mat V = arma::randg<arma::mat>(D,L);
  arma::mat logV = log(V);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  //arma::mat U = arma::zeros<arma::mat>(N,L);
  arma::vec z = arma::ones<arma::vec>(N);
  arma::vec ll(iter);
  double tau = 1.0;
  for (int i=0; i<iter; i++) {
    up_A(alpha, R, logV, y, xi, xp, a);
    up_B2(N, beta, V, logV, z, alpha, xi, xp, varind, b);
    up_tau(tau, z, y, V, xi,xp);
    ll.row(i) = lowerbound_logML2(z, alpha, beta, tau, V, logV, R, y, xi, xp, a, b);
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("precision")=tau,
                      Named("ELBO")=ll);
}

/*
List doVB_pois_missing(const arma::vec & y,
                       const arma::uvec & xi,
                       const arma::uvec & xp,
                       const arma::uvec & miss_row,
                       const arma::uvec & miss_col,
                       const arma::vec & sumx,
                       const arma::uvec & varind,
                       const int & D,
                       const int & L,
                       const int & iter,
                       const double & a,
                       const double & b,
                       const double & eta){
  int N = y.n_rows;
  arma::mat lambda = arma::randg<arma::mat>(D,L);
  arma::vec etahat = sumx+eta;
  arma::mat Xprob = arma::zeros<arma::mat>(miss_row.n_rows, D);
  arma::mat loglambda = log(lambda);
  arma::mat alpha = arma::ones<arma::mat>(D, L);
  arma::mat beta  = arma::ones<arma::mat>(D, L);
  arma::vec R = arma::zeros<arma::vec>(N);
  arma::vec logw = arma::ones<arma::vec>(D);
  logw /= sum(logw);
  logw = log(logw);
  arma::vec ll(iter);
  arma::vec ll2(iter);
  for (int i=0; i<iter; i++) {
    up_x(Xprob, y, xi, xp, sumx, miss_row, miss_col, varind, loglambda, lambda, D, L, logw, etahat, eta);
    arma::mat r =  myprod(N, xi, xp, exp(loglambda));
    r.rows(miss_row) %=  exp(Xprob*loglambda);
    R = sum(r, 1);
    alpha = mysum_t(D, xi, xp, r.each_col()%(y/R)) + a;
    up_B(N, beta, lambda, loglambda, alpha, xi, xp, varind, b);
    ll.row(i) = lowerbound_logML_pois(alpha, beta, lambda, loglambda, R, y, a, b);
    ll2.row(i) = lowerbound_logML_mult(logw, sumx, Xprob, etahat, eta);
  }
  return List::create(Named("shape")=alpha,
                      Named("rate")=beta,
                      Named("Xprob")=Xprob,
                      Named("shape_x")=etahat,
                      Named("ELBO")=ll, Named("ELBO2")=ll2);
}
*/