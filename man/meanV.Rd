\name{meanVy}
\alias{meanV}
\title{Compute Mean of Posterior Distributions}
\description{
  Computes the mean of Gamma-distributed latent variables.
}
\usage{
meanV(out,  logarithm = FALSE)
}
\arguments{
  \item{out}{List containing \code{shape} and \code{rate} matrices.}
  \item{logarithm}{If \code{TRUE}, return the mean of the logarithms, which is different from the logarithm of the mean.}
}
\examples{
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

V = meanV(out)
}