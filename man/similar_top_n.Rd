\name{similar_top_n}
\alias{similar_top_n}
\title{
get top n elements that have higher similarity from a matrix
}
\description{
get top \code{n} elements that have higher similarity about a vector \code{ref} from a matrix \code{V}
}
\usage{
simmilar_top_n(V, ref,
               n=10L,
               decreasing=FALSE)
}

\arguments{
  \item{V}{
target matrix
}
\item{ref}{
refference vector
}
\item{n}{
size of the return value
}
\item{decreasing}{
If \code{TRUE}, return top n elements which have higher dissimilarity
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
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
