\name{prod_mPPCA}
\alias{prod_mPPCA}
\alias{prod_mPPCA.default}
\alias{prod_mPPCA.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix product
}
\description{
\eqn{\prod_{d=1}^D\mu_{dl}^{x_{nd}}}
}
\usage{
prod_mPPCA.formula(formula,
                   data = parent.frame(), mu)


prod_mPPCA.default(X, mu)


prod_mPPCA(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  sparse design matrix
}
  \item{mu}{
  matrix of parameter
}
  \item{formula}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{	
  an optional data frame, containing the variables in the formula.
}
  \item{\dots}{
  additional arguments to be passed to the low level fitting functions 
}
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
