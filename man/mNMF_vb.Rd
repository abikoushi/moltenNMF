\name{mNMF_vb}
\alias{mNMF_vb}
\alias{mNMF_vb.default}
\alias{mNMA_vb.formula}
\title{
Estimate the model of molten NMF using variational Bayes
}
\description{
\code{mNMF_vb} is used to fit molten NMF.
}
\usage{
mNMF_vb.default(y,
                X,
                L,
                iter=1000, 
                a=0.5, 
                b=0.01,
                V=NULL,
                offset = NULL,
                display_progress=TRUE,
                indices=NULL)


mNMF_vb.formula(formula,
                data = parent.frame(),
                L,
                iter=1000,
                a=0.5,
                b=0.01,
                V=NULL,
                offset = NULL,
                display_progress=TRUE)


mNMF_vb(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  response variable (nonnegative-integer)
}
  \item{X}{
  sparse design matrix
}
  \item{formula}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{	
  an optional data frame, containing the variables in the formula.
}
  \item{L}{
  number of componets
}
  \item{iter}{
  number of iteration
}
  \item{a}{
  shape parameter of gamma prior distribution
}
  \item{b}{
  rate parameter of gamma prior distribution
}
  \item{V}{
  initial values of latent variable V. If \code{NULL} it choosed randomly. 
}
  \item{offset}{
  The vector of offset factor which scale output. \code{NULL} is equivalent to the vector of 1s.
  }
  \item{display_progress}{
  if \code{TRUE} show progress bar
}
 \item{indices}{
 The variables group indices. If design matrix \code{X} is made by \code{\link{sparse oneshot}}, keep default value.
 }
  \item{\dots}{
  additional arguments to be passed to the low level fitting functions
}
}
\value{
 \item{shape}{shape of the variational posterior (gamma distribution)}
 \item{rate}{rate of the variational posterior (gamma distribution)}
 \item{ELBO}{evidence lower bound}
}
\seealso{
\code{\link{mrNMF_vb}}
}
\examples{

Titanicdf <- as.data.frame(Titanic)
set.seed(794)
out2 <- mNMF_vb(Freq~Survived+Class+Sex+Age, data=Titanicdf, L=2, iter=500)

#Check convergence
plot(out2$ELBO, type = "l", xlab="iter", ylab="ELBO") 

#Check goodness-of-fit
V <- out2$shape/out2$rate
yhat <- product_m(~Survived+Class+Sex+Age, data=Titanicdf,V)
plot(log1p(yhat), log1p(Titanicdf$Freq), ylab = "fit", xlab = "obs", xlim = c(-0.5,7), ylim=c(-0.5,7))
abline(0, 1, lty=2, col="steelblue")

#see also: 
#https://github.com/abikoushi/moltenNMF/blob/main/Rmd/moltenNMF.md
}
