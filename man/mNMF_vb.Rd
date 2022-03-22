\name{mNMF_vb}
\alias{mNMF_vb}
\alias{mNMF_vb.default}
\alias{mNMA_vb.formula}
\title{
Estimates the model of molten NMF using variational Bayes
}
\description{
\code{mNMF_vb} is used to fit molten NMF.
}
\usage{
mNMF_vb.default(y,
                X,
                L,
                iter=1000, 
                a=1, 
                b=1)


mNMF_vb.formula(formula,
                data = parent.frame(),
                L,
                iter=1000,
                a=1,
                b=1)


mNMF_vb(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  response variable
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
  \item{X}{
  number of componets
}
  \item{iter}{
  number of iteration
}
  \item{a}{
  shape parameter of gamma prior distribution
}
  \item{a}{
  rate parameter of gamma prior distribution
}
  \item{\dots}{
  additional arguments to be passed to the low level fitting functions 
}
}
\value{
 \item{shape }{shape of the variational posterior (gamma distribution)}
 \item{rate }{rate of the variational posterior (gamma distribution)}
 \item{ELBO }{evidence lower bound}
 \item{indices }{indices}
}
\author{
Ko ABE
}
\seealso{
\code{\link{mPPCA_vb}}
}
\examples{
### see vignette
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
