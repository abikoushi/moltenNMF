\name{rowmeanvar_mtx}
\alias{rowmeanvar_mtx}
\title{Compute Row-wise Mean and Variance from Sparse Text Matrix}
\description{
Computes the mean and variance of each row from a sparse matrix stored in a text file.
}
\usage{
rowmeanvar_mtx(file_path, n_header=2)
}
\arguments{
  \item{file_path}{Path to the text file containing non-zero entries in "i j v" format.}
  \item{n_header}{Number of header lines to skip.}
}
\details{
The function reads a sparse matrix in coordinate (COO) format and accumulates values for each row to compute sample mean and variance.  
It assumes that missing values are implicitly zero.

Variance is computed as:
\deqn{\mathrm{Var}(X_i) = \frac{1}{n-1} \sum_j (X_{ij} - \bar{X}_i)^2}
}
\value{
A list containing:
  \item{mean}{A numeric vector of row means.}
  \item{var}{A numeric vector of row variances.}
}
\examples{
## Not run:
library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "obs.mtx")
rowmeanvar_mtx("obs.mtx")
## End(**Not run**)
}
