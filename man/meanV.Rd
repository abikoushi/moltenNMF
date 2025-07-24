\name{meanVy}
\alias{meanV}
\title{Compute Mean of Posterior Distributions}
\description{
  Computes the mean of Gamma-distributed latent variables.
}
\usage{
meanV(out,  logarithm = FALSE)
}
\arguments{
  \item{out}{List containing \code{shape} and \code{rate} matrices.}
  \item{logarithm}{If \code{TRUE}, return the mean of the logarithms, which is different from the logarithm of the mean.}
}
\examples{
Titanicdf <- as.data.frame(Titanic)
set.seed(794)
out2 <- mNMF_vb(Freq~Survived+Class+Sex+Age, data=Titanicdf, L=2, iter=500)
head( meanV(out2) ) #equivalent to follows:
head(out2$shape/out2$rate)
}
