\name{rpredictor_mNMF}
\alias{dpredictor_mNMF}
\alias{rpredictor_mNMF}
\alias{rpredictor_mrNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampler for predictive distribution of mNMF.
}
\description{
Get sample frrom posterior predictive distribution from the fitted model of mNMF using random number.
}
\usage{
rpredictor_mNMF(X, np, shape, rate)

rpredictor_mrNMF(X, np, shape, rate, precision)

dpredictor_mNMF(y, X, np, shape, rate)

dpredictor_mrNMF(y, X, np, shape, rate, precision)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  response variable (nonnegative-integer)
}
  \item{X}{
  sparse design matrix
}
  \item{np}{
  size of the sample (Natural number)
}
  \item{shape}{
  shape of the variational posterior (gamma distribution)
}
  \item{rate}{rate of the variational posterior (gamma distribution)}
 
 \item{precision}{the parameter of prior distribution (gamma distribution)}
}

\seealso{
\code{\link{sparse_onehot}}
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
