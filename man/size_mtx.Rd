\name{size_mtx}
\alias{size_mtx}
\title{Read Matrix Dimensions from File}
\description{
  Reads matrix dimensions from a file in Matrix Market format or similar structure.
}
\usage{
size_mtx(file_path)
}
\arguments{
  \item{file_path}{Path to the input file containing matrix data.}
}
\details{
  This function reads the second line of the file (after skipping the Matrix Market header), assuming it contains the matrix dimensions.
}
\value{
  An integer vector containing the matrix dimensions; the number of rows, columns, and nonzero elements.
}
\examples{
## Not run:
library(Matrix)
i <- c(1,3:8); j <- c(2,9,6:10); x <- 7 * (1:7)
(A <- sparseMatrix(i, j, x = x))  ##  8 x 10 "dgCMatrix"
writeMM(obj = A, file = "obs.mtx")
size_mtx("obs.mtx") #[1]  8 10  7
## End(**Not run**)
}
