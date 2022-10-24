\name{sparse_onehot}
\alias{sparse_onehot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
make sparse one-hot encoding (design matrix)
}
\description{
make sparse one-hot encoding (design matrix)
}
\usage{
sparse_onehot(object,
              data = environment(object),
              xlev = NULL,
              sep = "_",
              na.action='na.pass',
              to = "l")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{
  an optional data frame, containing the variables in the formula.
}
}
\details{
make sparse one-hot encoding (design matrix)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##see vignette
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
