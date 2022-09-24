\name{product_m}
\alias{product_m}
\alias{product_m.default}
\alias{product_m.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix product
}
\description{
calculate \eqn{f_{nl} = \prod_{d=1}^Dv_{dl}^{x_{nd}}}. If `aggregate = TRUE`, return \eqn{\sum_{l=1}^L f_{nl}}.
}
\usage{
product_m.formula(formula,
  data = parent.frame(),
  V,
  aggregate = TRUE)


product_m.default(X,
  V,
  aggregate = TRUE)


product_mNMF(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  sparse design matrix
}
  \item{V}{
  matrix of parameter
}
  \item{formula}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{	
  an optional data frame, containing the variables in the formula.
}
  \item{aggregate}{	
  If `TRUE`, return sum of the values of components. 
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
