\name{NMF2D_svb_bin}
\alias{NMF2D_svb_bin}
\title{Stochastic Variational Bayes for Poisson NMF with Binary Input}
\description{
  Stochastic VB algorithm for large-scale matrix factorization using binary files.
}
\usage{
NMF2D_svb_bin(filepath_x, filepath_y, dims, rank, 
              n_epochs, n_baches, lr_param, 
              lr_type = "exponential", 
              prior_shape = 1, prior_rate = 1, 
              display_progress = TRUE)
}
\arguments{
  \item{filepath_x}{Path to binary file for row and column indices.}
  \item{filepath_y}{Path to binary file for observed values.}
  \item{dims}{Vector of dimensions: (number of rows, columns, and non-zero observations).}
  \item{rank}{Number of latent components.}
  \item{n_epochs}{Number of training epochs.}
  \item{n_baches}{Number of mini-batches.}
  \item{lr_param}{Learning rate parameters. For \code{exponential},
      \deqn{\rho = (\mathtt{lr\_param[1]} + t)^{-\mathtt{lr\_param[2]}}}{
    rho = (lr_param[1] + t)^(-lr_param[2])
    }
    where \eqn{t} is the current iteration number.
    For \code{const}, it should be a constant learning rate value.
    }
  \item{lr_type}{A character string indicating the learning rate schedule.
  Only \code{"exponential"} and \code{"const"} are currently implemented.}
  \item{prior_shape}{Shape parameter for Gamma prior.}
  \item{prior_rate}{Rate parameter for Gamma prior.}
  \item{display_progress}{Logical; whether to display progress.}
}
\value{
  Posterior parameters for the factor matrices.
}
\seealso{
  \code{\link{NMF2D_vb}}, \code{\link{NMF2D_svb_bin}}
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
moltenNMF:::writebin_spmat(dat$Y,
                           filepath_x = "test_x.bin", 
                           filepath_y = "test_y.bin")
dims = c(dim(dat$Y), nnzero(dat$Y))            
system.time({
  res = moltenNMF:::NMF2D_svb_bin(filepath_x = "test_x.bin",
                                  filepath_y = "test_y.bin",
                                  dims = dims,
                                  rank=3, n_epochs=200, n_baches=1000,
                                  lr_param = c(16, 0.8),
                                  lr_type = "exponential",
                                  display_progress = TRUE)  
})

V = moltenNMF:::meanV_array(res)
V = moltenNMF:::rearrange_cols(V, normalize = FALSE)
fit2 = V[[1]]\%*\%t(V[[2]])
plot(x=c(fit2), y=c(as.matrix(dat$Y)), pch=".", xlab="fitted", ylab="observed")
abline(0,1,col="lightgrey")
## End(**Not run**)
}