\name{obsfitloss_2d_mtx}
\alias{obsfitloss_2d_mtx}
\title{Compute Poisson and MSE Loss from Two Factor Matrices}
\description{
Calculates Poisson loss and Mean Squared Error (MSE) between observed sparse entries and reconstructed values from two low-rank factor matrices.
}
\usage{
obsfitloss_2d_mtx(readtxt, V1, V2, n_header)
}
\arguments{
  \item{readtxt}{Path to a sparse matrix text file. Each line contains space-separated integers: row index, column index, and value.}
  \item{V1}{A numeric matrix representing the first factor matrix (e.g., rows × rank).}
  \item{V2}{A numeric matrix representing the second factor matrix (e.g., columns × rank).}
  \item{n_header}{Number of header lines to skip in the text file.}
}
\details{
Each non-zero observation from the text file is compared to the reconstructed value computed as the inner product of the corresponding rows in \code{V1} and \code{V2}.  
The losses are averaged over all possible entries (i.e., full matrix of size \code{nrow(V1)} × \code{nrow(V2)}).

Poisson loss is computed as:
\deqn{\frac{1}{N} \sum (y_{ij} \cdot \log \hat{y}_{ij} - \hat{y}_{ij})}

MSE is computed as:
\deqn{\frac{1}{N} \sum (y_{ij} - \hat{y}_{ij})^2}
}
\value{
A list with:
  \item{Poisson}{Mean Poisson loss.}
  \item{MSE}{Mean Squared Error.}
}
\examples{
V1 <- matrix(rexp(50), 10, 5)
V2 <- matrix(rexp(60), 12, 5)
loss <- obsfitloss_2d_mtx("obs.txt", V1, V2, n_header = 2)
print(loss)
}
\seealso{
\code{\link{obsfitloss_mtx}}
}
