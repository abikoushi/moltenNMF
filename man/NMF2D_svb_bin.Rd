\name{NMF2D_svb_bin}
\alias{NMF2D_svb_bin}
\title{Stochastic Variational Bayes for Poisson NMF with Binary Input}
\description{
  Stochastic VB algorithm for large-scale matrix factorization using binary files.
}
\usage{
NMF2D_svb_bin(filepath_x, filepath_y, dims, rank, 
              n_epochs, n_baches, lr_param, 
              lr_type = "exponential", subiter = 1, 
              prior_shape = 1, prior_rate = 1, 
              Vini = NULL, display_progress = TRUE)
}
\arguments{
  \item{filepath_x}{Path to binary file for row and column indices.}
  \item{filepath_y}{Path to binary file for observed values.}
  \item{dims}{Vector of dimensions: (number of rows, columns, and non-zero observations).}
  \item{rank}{Number of latent components.}
  \item{n_epochs}{Number of training epochs.}
  \item{n_baches}{Number of mini-batches.}
  \item{lr_param}{Learning rate parameters.}
  \item{lr_type}{Type of learning rate schedule.}
  \item{subiter}{Sub-iteration count (not used currently).}
  \item{prior_shape}{Shape parameter for Gamma prior.}
  \item{prior_rate}{Rate parameter for Gamma prior.}
  \item{Vini}{Initial values (not currently used).}
  \item{display_progress}{Logical; whether to display progress.}
}
\value{
  Posterior parameters for the factor matrices.
}
