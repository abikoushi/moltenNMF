\name{sparse_model_matrix_b}
\alias{sparse_model_matrix_b}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create sparse design matrix
}
\description{
This is a wrapper function of \code{sparse.model.matrix}.
}
\usage{
sparse_model_matrix_b(object, data = environment(object), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{
  an optional data frame, containing the variables in the formula.
}
  \item{\dots}{
  additional arguments to be passed to \code{sparse.model.matrix}
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
