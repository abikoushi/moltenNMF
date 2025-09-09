\name{rowfilter_mtx}
\alias{rowfilter_mtx}
\title{Filter Rows from a Sparse Matrix File and Write to New File}
\description{
Reads a sparse matrix text file and writes a filtered version containing only specified rows to a new file.
}
\usage{
rowfilter_mtx(readtxt,
              writetxt,
              rowind,
              n_header = 2)
}
\arguments{
  \item{readtxt}{Path to the input sparse matrix text file.}
  \item{writetxt}{Path to the output file to be created.}
  \item{rowind}{An integer vector of row indices (1-based) to retain.}
  \item{n_header}{Number of header lines in \code{readtxt}.}
}
\details{
This function preserves only the rows specified in \code{rowind}, adjusts the row indices to match the new index (starting from 1), and updates the matrix size in the header line accordingly.

Assumes the text file is in "i j v" format with a \code{n_header}-line header:
\itemize{
  \item Line 1: comment or metadata (copied unchanged)
  \item Line 2: original matrix size and non-zero count (replaced in output)
}
}
\value{
No return value. A new file is written to \code{writetxt}.
}
\examples{
## Not run:
library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "obs.mtx")
# Keep only rows 1, 5
rowfilter_mtx("obs.mtx", "filtered.mtx", c(1, 5))
readMM("filtered.mtx")
## End(**Not run**)
}
