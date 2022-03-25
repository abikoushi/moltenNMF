\name{mBMF_mcvb}
\alias{mBMF_mcvb}
\alias{mBMF_mcvb.default}
\alias{mBMF_mcvb.formula}
\title{
Estimate the model of molten binary matrix factorization using Monte-Carlo variation Bayes
}
\description{
\code{mBMF_mcvb} is used to fit molten binary matrix factorization.
}
\usage{
mBMF_mcvb.default(y, 
                   X,
                   L,
                   iter=2000,
                   lambda=1,
                   tau=1)
                    
                    
mBMF_mcvb.formula(formula,
                   data = parent.frame(),
                   L,
                   iter=2000,
                   lambda=1,
                   tau=1)


mBMF_mcvb(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  response variable (0 or 1)
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
  \item{lambda}{
  precision parameter (inverse of variance)
}
  \item{tau}{
  precision parameter of prior distribution
}
  \item{\dots}{
  additional arguments to be passed to the low level fitting functions 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{mean }{mean of the variational posterior (normal distribution)}
 \item{sd }{standard deviation of the variational posterior (normal distribution)}
 \item{logloss }{log-probabilty evalueted by current model}
 \item{indices }{indices of variavles}
}
\author{
ko ABE
}
\seealso{
\code{\link{mNMF_vb}}
}
\examples{
###
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
