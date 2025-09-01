\name{rearrange_cols}
\alias{rearrange_cols}
\title{Reorder Latent Components}
\description{
  Reorders columns of factor matrices based on variability or other summary statistics.
}
\usage{
rearrange_cols(Vm, axis = 1L, FUN = var, normalize = FALSE, decreasing = TRUE)
}
\arguments{
  \item{Vm}{List of factor matrices.}
  \item{axis}{Which matrix to use for determining the order.}
  \item{FUN}{Function used to compute importance (e.g., \code{var}, \code{mean}).}
  \item{normalize}{Logical; whether to normalize each row.}
  \item{decreasing}{Logical; whether to sort in decreasing order.}
}
\value{
  List of reordered factor matrices.
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
set.seed(1234)
dat <- set_data_mf(3, 99, 200)
system.time({
  out <- moltenNMF:::NMF2D_vb(dat$Y, rank = 3, iter = 1000)
})
V = meanV(out)
head(V[[1]])
V = rearrange_cols(V)
head(V[[1]])
## End(**Not run**)
}