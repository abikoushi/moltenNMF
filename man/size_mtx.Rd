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
  An integer vector containing the matrix dimensions.
}
\seealso{
  \code{\link{scan}}, \code{\link{file}}
}
