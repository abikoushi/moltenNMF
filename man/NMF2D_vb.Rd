\name{NMF2D_vb}
\alias{NMF2D_vb}
\title{Variational Bayes for Poisson NMF (for 2D array format)}
\description{
  Applies Variational Bayes (VB) inference for Poisson non-negative matrix factorization (NMF) on sparse matrix.
}
\usage{
NMF2D_vb(Y, rank, iter = 100, weight = NULL, 
         prior_shape = 1, prior_rate = 1, 
         Vini = NULL, dims = NULL, 
         display_progress = TRUE)
         
         
}
\arguments{
  \item{Y}{A sparse matrix (class \code{dgTMatrix} or coercible to it).}
  \item{rank}{Integer, the number of latent components.}
  \item{iter}{Number of VB iterations.}
  \item{weight}{Optional list of row/column weights.}
  \item{prior_shape}{Shape parameter of the Gamma prior.}
  \item{prior_rate}{Rate parameter of the Gamma prior.}
  \item{Vini}{Optional list of initial factor matrices.}
  \item{dims}{Optional vector specifying matrix dimensions.}
  \item{display_progress}{Logical; whether to print progress.}
}
\value{
  A list containing posterior shape and rate parameters for each factor matrix.
}
\seealso{
  \code{\link{NMF2D_svb}}, \code{\link{meanV_array}}, \code{\link{rearrange_cols}}
}
\examples{
library(Matrix)
Y <- rsparsematrix(100, 80, density = 0.1)
res <- NMF2D_vb(Y, rank = 5)
}
