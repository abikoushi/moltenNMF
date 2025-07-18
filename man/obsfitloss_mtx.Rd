\name{obsfitloss_mtx}
\alias{obsfitloss_mtx}
\title{Compute Poisson and MSE Loss from Observed Sparse Matrix File}
\description{
Calculates the Poisson loss and Mean Squared Error (MSE) between observed sparse data and a predicted dense matrix.
}
\usage{
obsfitloss_mtx(readtxt, fit, n_header)
}
\arguments{
  \item{readtxt}{A string giving the file path to a sparse matrix text file. Each line should contain three space-separated values: row index, column index, and value (i j v).}
  \item{fit}{A numeric dense matrix representing predicted values.}
  \item{n_header}{An integer specifying the number of header lines to skip in the text file.}
}
\details{
This function reads a sparse matrix from a text file and computes losses between the observed entries and the corresponding values in the prediction matrix \code{fit}. Observed entries are zeroed out in \code{fit} before computing loss over the remaining unobserved entries.

Poisson loss is computed as:
\deqn{- (x \log y - y - \log(x!))}

Mean Squared Error (MSE) is computed over all elements (observed and unobserved).
}
\value{
A list with the following elements:
  \item{Poisson}{Mean Poisson loss over all entries.}
  \item{MSE}{Mean Squared Error over all entries.}
}
\examples{
## Not run:
fit <- matrix(runif(100), nrow=10)
loss <- obsfitloss_mtx("obs.mtx", fit, n_header=2)
print(loss)
## End(**Not run**)
}
\seealso{
\code{\link{obsfitloss_2d_mtx}}
}
