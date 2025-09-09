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
  \code{\link{NMF2D_svb}}
}
\examples{
## Not run:
library(Matrix)
set_data_mf <- function(L, nrow, ncol, mu=0){
  W <- matrix(rnorm(nrow*L,0,1), ncol=L)
  H <- matrix(rnorm(L*ncol,0,1), nrow=L)
  W <- sweep(W,1,rowMeans(W)-mu)
  H <- sweep(H,2,colMeans(H)-mu)
  Y <- matrix(rpois(nrow*ncol, exp(W)\%*\%exp(H)), nrow, ncol)
  Y <- as(Y, "TsparseMatrix")
  list(Y=Y, trueW=W, trueH=H)
}
dat <- set_data_mf(3, 99, 500)
system.time({
  out <- moltenNMF:::NMF2D_vb(dat$Y, rank = 3, iter = 1000)
})
plot(out$ELBO[-1], type = "l")
## End(**Not run**)
}