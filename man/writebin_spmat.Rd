\name{writebin_spmat}
\alias{writebin_spmat}
\title{Write Sparse Matrix to Binary Files}
\description{
  Writes a sparse matrix in triplet form (\code{dgTMatrix}) to binary files representing indices and values.
}
\usage{
writebin_spmat(object, filepath_x, filepath_y, index_decrement = 0L)
}
\arguments{
  \item{object}{A sparse matrix of class \code{"dgTMatrix"} or coercible to it.}
  \item{filepath_x}{File path to write the concatenated row and column indices.}
  \item{filepath_y}{File path to write the non-zero values.}
  \item{index_decrement}{An integer to subtract from each index (default is \code{0L}). Useful for converting to 0-based indexing.}
}
\details{
  The function writes the non-zero values to \code{filepath_y} using \code{writeBin}, and the row/column indices (interleaved) to \code{filepath_x}. If the matrix is not in triplet format, it is converted.
}
\seealso{
  \code{\link{Matrix}}, \code{\link{writeBin}}
}
