#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
List rowmeanvar_mtx(const int & n_row, const int & n_col,
                    const std::string & readtxt,
                    const int & n_header) {
  int x;
  int y;
  double v;
  std::ifstream file(readtxt);
  std::string str;    
  int index = 0;
  int n = 0;
  arma::vec vout = arma::zeros<arma::vec>(n_row);
  arma::vec v2out = arma::zeros<arma::vec>(n_row);
  for(int i=0; i<n_header; i++){
    std::getline(file, str);
    index++;
  }
  //Rprintf("%d\n", index);
  while(std::getline(file, str)){
    //if(index > 2){
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
    vout(x-1) += v;
    v2out(x-1) += pow(v,2);
    //}
    index++;
  }
  vout /= n_col;
  v2out /= n_col-1;
  v2out -= pow(vout,2);
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
