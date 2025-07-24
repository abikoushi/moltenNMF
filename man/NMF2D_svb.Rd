\name{NMF2D_svb}
\alias{NMF2D_svb}
\title{Stochastic Variational Bayes for Poisson NMF (for 2D array format)}
\description{
  Applies stochastic variational inference for Poisson NMF using stochastic mini-batches.
}
\usage{
NMF2D_svb(Y, rank, n_epochs, n_baches, 
          lr_param = c(1, 0.8), lr_type = "exponential",
          dims = NULL, prior_shape = 1, prior_rate = 1, 
          Vini = NULL, weight = NULL, 
          display_progress = TRUE)
}
\arguments{
  \item{Y}{A sparse matrix (class \code{dgTMatrix} or coercible to it).}
  \item{rank}{Number of latent components.}
  \item{n_epochs}{Number of passes over the data.}
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
  \item{dims}{Optional vector specifying matrix dimensions.}
  \item{prior_shape}{Shape of the Gamma prior.}
  \item{prior_rate}{Rate of the Gamma prior.}
  \item{Vini}{Initial factor matrices.}
  \item{weight}{Optional list of row/column weights.}
  \item{display_progress}{Logical; whether to show progress bar.}
}
\value{
  A list with variational parameters learned through stochastic inference.
}
\seealso{
  \code{\link{NMF2D_vb}}, \code{\link{NMF2D_svb_bin}}
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
  out2 <- moltenNMF:::NMF2D_svb(dat$Y, rank = 3,
                                n_epochs = 200, n_baches = as.integer(2000),
                                prior_shape = 1, prior_rate = 1,
                                lr_param = c(15, 0.8),
                                lr_type = "exponential")
})
plot(out2$ELBO[-1], type="l")
}
