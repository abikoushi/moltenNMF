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
###
}