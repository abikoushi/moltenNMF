\name{slice_rows}
\alias{slice_rows}
\title{
Slice Rows from a Design Matrix
}
\description{
Extracts specific rows from a design matrix while preserving relevant attributes.
}
\usage{
slice_rows(x, i)
}
\arguments{
  \item{x}{A design matrix which rows will be extracted.}
  \item{i}{An index vector specifying the rows to be selected.}
}
\details{
This function slices the specified rows from a given matrix or data frame, ensuring that
certain attributes such as \code{indices}, \code{term.labels}, \code{value.labels}, and \code{assign}
are preserved in the output when they exist in the input.
}

\seealso{
\code{\link{sparse_onehot}}
}
\examples{
df <- expand.grid(x1=factor(1:10), x2=factor(1:10))
X <- sparse_onehot(~x1+x2, data=df)
sub_x <- slice_rows(X, c(1, 3, 5))
}
