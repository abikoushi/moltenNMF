\name{rowmeanvar_mtx}
\alias{rowmeanvar_mtx}
\title{Compute Row-wise Mean and Variance from Sparse Text Matrix}
\description{
Computes the mean and variance of each row from a sparse matrix stored in a text file.
}
\usage{
rowmeanvar_mtx(n_row, n_col, readtxt, n_header)
}
\arguments{
  \item{n_row}{Total number of rows in the full matrix.}
  \item{n_col}{Total number of columns in the full matrix.}
  \item{readtxt}{Path to the text file containing non-zero entries in "i j v" format.}
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
rowmeanvar_mtx(n_row = 100, n_col = 500, readtxt = "obs.mtx", n_header = 2)
}
\seealso{
\code{\link{rowfilter_mtx}}
}
