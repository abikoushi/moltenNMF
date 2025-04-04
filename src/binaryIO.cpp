#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "binaryIO.h"
using namespace Rcpp;

//////
//read
//rows: rows to read
//Y, X:template matrix
//////
void readRowsFromBinary(arma::mat & Y, 
                        const std::string & filepath,
                        const arma::uvec & rows) {
  std::ifstream in(filepath, std::ios::binary);
  if (!in.is_open()) {
    Rcpp::stop("Failed to open file for reading.");
  }
  int nrows = rows.n_elem;
  for (arma::uword i = 0; i < nrows; ++i) {
    //Rprintf("%d ", i);
    size_t offset = rows(i) * Y.n_cols * sizeof(double);
    in.seekg(offset, std::ios::beg);
    if (!in) {
      Rcpp::stop("Error seeking to the specified position in the file.");
    }
    std::vector<double> buffer(Y.n_cols);
    in.read(reinterpret_cast<char*>(buffer.data()), Y.n_cols * sizeof(double));
    if (!in) {
      Rcpp::stop("Error reading from the file.");
    }
    Y.row(i) = arma::rowvec(buffer);
  }
  in.close();
}

void readRowsFromBinary_umat(arma::umat & X, 
                             const std::string & filepath,
                             const arma::uvec & rows) {
  std::ifstream in(filepath, std::ios::binary);
  if (!in.is_open()) {
    Rcpp::stop("Failed to open file for reading.");
  }
  int nrows = rows.n_elem;
  for (arma::uword i = 0; i < nrows; ++i) {
    size_t offset = rows(i) * X.n_cols * sizeof(unsigned int);
    in.seekg(offset, std::ios::beg);
    if (!in) {
      Rcpp::stop("Error seeking to the specified position in the file.");
    }
    std::vector<unsigned int> buffer(X.n_cols);
    in.read(reinterpret_cast<char*>(buffer.data()), X.n_cols * sizeof(unsigned int));
    if (!in) {
      Rcpp::stop("Error reading from the file.");
    }
    for (arma::uword j = 0; j < X.n_cols; ++j) {
      X(i, j) = buffer[j];
    }
  }
  in.close();
}

// [[Rcpp::export]]
List read_bin(const std::string & filepath_x,
              const std::string & filepath_y,
              const arma::uvec & bag,
              const int & x_dim){
  arma::umat X(bag.n_rows, x_dim);
  arma::vec val(bag.n_rows);
  readRowsFromBinary_umat(X, filepath_x, bag);
  readRowsFromBinary(val, filepath_y, bag);
  return List::create(_["X"]=X, _["y"]=val);
}

///
//write
///

/*
void writeBinaryVec(arma::vec data, std::string filepath) {
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    Rcpp::stop("Failed to open.");
  }
  out.write(reinterpret_cast<char*>(data.memptr()), data.n_elem * sizeof(double));
  out.close();
}

void writeBinaryFile_umat(arma::umat data, std::string filepath) {
  std::ofstream out(filepath, std::ios::binary);
  if (!out.is_open()) {
    Rcpp::stop("Failed to open.");
  }
  out.write(reinterpret_cast<char*>(data.memptr()), data.n_elem * sizeof(unsigned int));
  out.close();
}
*/