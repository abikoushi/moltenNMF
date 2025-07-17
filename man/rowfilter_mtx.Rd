\name{rowfilter_mtx}
\alias{rowfilter_mtx}
\title{Filter Rows from a Sparse Matrix File and Write to New File}
\description{
Reads a sparse matrix text file and writes a filtered version containing only specified rows to a new file.
}
\usage{
rowfilter_mtx(readtxt, writetxt, rowind)
}
\arguments{
  \item{readtxt}{Path to the input sparse matrix text file.}
  \item{writetxt}{Path to the output file to be created.}
  \item{rowind}{An integer vector of row indices (1-based) to retain.}
}
\details{
This function preserves only the rows specified in \code{rowind}, adjusts the row indices to match the new index (starting from 1), and updates the matrix size in the header line accordingly.

Assumes the text file is in "i j v" format with a 2-line header:
\itemize{
  \item Line 1: comment or metadata (copied unchanged)
  \item Line 2: original matrix size and non-zero count (replaced in output)
}
}
\value{
No return value. A new file is written to \code{writetxt}.
}
\examples{
# Keep only rows 1, 5, and 10
rowfilter_mtx("obs.mtx", "filtered.txt", c(1, 5, 10))
}
