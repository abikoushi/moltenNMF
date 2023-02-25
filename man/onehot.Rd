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
              na.action='na.pass')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{
  an optional data frame, containing the variables in the formula.
}
\item{xlev}{
a named list of character vectors giving the full set of levels to be assumed for each factor.
}
\item{sep}{
character of separator within \code{colnames}
}
\item{na.action}{
how \code{NA}s are treated (see \code{\link{model.frame}})
}
}
\details{
make sparse one-hot encoding (design matrix)
}

\examples{
set.seed(1234)
df <- data.frame(x1 = sample(letters[1:3], 10, replace = TRUE),
                 x2 = sample(LETTERS[1:3], 10, replace = TRUE))
sparse_onehot(~x1+x2, data=df)

sparse_onehot(~x1+x2, data=df, dummy_m=2)
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
