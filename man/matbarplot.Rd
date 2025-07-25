\name{matbarplot}
\alias{matbarplot}

\title{
Simple barplot for nonnegative matrix
}
\description{
This is a function that easily illustrates a non-negative matrix.
To see the data in more detail, we recommend to try another package.
}
\usage{
matbarplot(m,
  xgroup=NULL,
  col=hcl.colors,
  ...)
}
\arguments{
  \item{m}{
matrix for plot
}
  \item{xgroup}{
annotation for x-axis
}
  \item{col}{
color palette (see \code{\link{hcl.colors}})
}
  \item{\dots}{
other arguments for barplot function
}
}
\examples{
matbarplot(as.matrix(iris[,-5]), xgroup = iris$Species)
}
