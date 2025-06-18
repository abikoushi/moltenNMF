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
