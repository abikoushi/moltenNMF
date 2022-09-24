\name{rpredictor_mNMF}
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
rpredictor_mNMF(X, np, alpha, beta)

rpredictor_mrNMF(X, np, alpha, beta, tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  sparse design matrix
}
  \item{np}{
  size of the sample (Natural number)
}
  \item{alpha}{
  shape of the variational posterior (gamma distribution)
}
  \item{beta}{rate of the variational posterior (gamma distribution)}
 
 \item{tau}{the parameter of prior distribution (gamma distribution)}
}

\seealso{
\code{\link{sparse_model_matrix_b}}
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
