#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include "logexpfuns.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace Rcpp;

double poisloss(const double & obs , const double & fit){
  return -(xlogy(obs, fit) - fit - lgamma(obs+1));
}

double mseloss(const double & obs , const double & fit){
  return pow(obs - fit, 2);
}

////
//post-
////
// [[Rcpp::export]]
List obsfitloss_mtx(const std::string & readtxt, arma::mat fit, const int & n_header){
  double MSE = 0;
  double pois = 0;
  int x;
  int y;
  double v;
  std::ifstream file(readtxt);
  std::string str;    
  int index = 0;
  // int n = 0;
  for(int i=0; i<n_header; i++){
    //skip header
    std::getline(file, str);
    index++;
  }
  while(std::getline(file, str)){
    //index++;
    std::stringstream ss(str);
    std::vector<std::string> svec;
    while( ss.good() ){
      std::string substr;
      getline(ss, substr, ' ');
      svec.push_back(substr);
    }
    // Rprintf("%d ", index);
    x = stoi(svec[0]);
    y = stoi(svec[1]);
    v = stod(svec[2]);
    x--;
    y--;
    pois += poisloss(v, fit(x,y));
    MSE += mseloss(v, fit(x,y));
    fit(x,y) = 0.;
  }
  MSE += accu(pow(fit,2));
  MSE /= (double) fit.n_elem;
  pois += accu(fit);
  pois /= (double) fit.n_elem;
  return List::create(Named("Poisson")=pois,
                      Named("MSE")=MSE);
}

// [[Rcpp::export]]
List obsfitloss_2d_mtx(const std::string & readtxt,
                       const arma::mat & V1,
                       const arma::mat & V2,
                       const int & n_header){
  double of = 0.0;
  double ologf = 0.0;
  double o2 = 0.0;
  int x;
  int y;
  double v;
  std::ifstream file(readtxt);
  std::string str;
  int index = 0;
  // int n = 0;
  for(int i=0; i<n_header; i++){
    //skip header
    std::getline(file, str);
    index++;
  }
  while(std::getline(file, str)){
    //index++;
    std::stringstream ss(str);
    std::vector<std::string> svec;
    while( ss.good() ){
      std::string substr;
      getline(ss, substr, ' ');
      svec.push_back(substr);
    }
    x = stoi(svec[0]);
    y = stoi(svec[1]);
    v = stod(svec[2]);
    x--;
    y--;
    double fit = arma::dot(V1.row(x), V2.row(y));
    of += v*fit;
    ologf += v*log(fit);
    o2 += pow(v,2);
  }
  double f2=0.0;
  double sumf =0.0;
  for(arma::uword i=0; i < V1.n_rows; i++){
    for(arma::uword j=0; j < V2.n_rows; j++){
      double fit = arma::dot(V1.row(i), V2.row(j));
      f2 += pow(fit,2);
      sumf += fit;
    } 
  }
  double N = V1.n_rows * V2.n_rows;
  double MSE = (o2+f2-2*of)/N;
  double pois =  (sumf - ologf)/N;
  return List::create(Named("Poisson")=pois,
                      Named("MSE")=MSE);
}

 
////
//pre-
////

// [[Rcpp::export]]
List rowmeanvar_txt(const int & n_row,
                    const int & n_col,
                    const std::string & readtxt,
                    const int & n_header) {
  int x;
  // int y;
  double v;
  std::ifstream file(readtxt);
  std::string str;
  int index = 0;
  // int n = 0;
  arma::vec vout = arma::zeros<arma::vec>(n_row);
  arma::vec v2out = arma::zeros<arma::vec>(n_row);
  for(int i=0; i<n_header; i++){
    //skip header
    std::getline(file, str);
    index++;
  }
  while(std::getline(file, str)){
    std::stringstream ss(str);
    std::vector<std::string> svec;
    while( ss.good() ){
      std::string substr;
      getline(ss, substr, ' ');
      svec.push_back(substr);
    }
    x = stoi(svec[0]);
    // y = stoi(svec[1]);
    v = stod(svec[2]);
    vout(x-1) += v;
    v2out(x-1) += pow(v,2);
    //index++;
  }
  vout /= (double) n_col;
  v2out /= ((double) n_col) - 1;
  v2out -= pow(vout, 2);
  return List::create(Named("mean")=vout,
                      Named("var")=v2out);
}

// [[Rcpp::export]]
void rowfilter_mtx(const std::string & readtxt,
                   const std::string & writetxt,
                   const arma::vec & rowind){
  int x;
  std::ifstream file(readtxt);
  std::string str;
  std::ofstream newfile;
  newfile.open(writetxt);
  // first line
  std::getline(file, str);
  newfile << str + "\n";
  // second line
  std::getline(file, str);
  std::stringstream ss(str);
  std::vector<std::string> svec;
  while( ss.good() ){
    std::string substr;
    getline(ss, substr, ' ');
    svec.push_back(substr);
  }
  std::string s0 = std::to_string(rowind.n_rows);
  s0 =  s0 + " " + svec[1];
  //myfile << "#size: "+ s0 + " " + svec[1] + " " + "\n";
  //filtering
  int nonzero = 0;
  std::vector<std::string> out;
  while (std::getline(file, str)){
      std::stringstream ss(str);
      std::vector<std::string> svec;
      while( ss.good() ){
        std::string substr;
        getline(ss, substr, ' ');
        svec.push_back(substr);
      }
      x = stoi(svec[0]);
      if(arma::any(rowind==x)){
        nonzero++;
        arma::uvec newx = find(rowind==x);
        std::string s0 = std::to_string(newx(0)+1);
        out.push_back(s0 + " " + svec[1] + " " + svec[2] + "\n");
        //myfile << s0 + " " + svec[1] + " " + svec[2] + "\n";
      }
  }
  //writing
  newfile << s0 + " " + std::to_string(nonzero) + "\n";
  for(int i=0; i<nonzero; i++){
    newfile << out[i];
  }
  newfile.close(); 
}
