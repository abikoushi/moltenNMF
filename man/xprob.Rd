\name{xprob_mNMF}
\alias{xprob_mNMF}
\alias{xprob_mNMF.formula}
\alias{xprob_mNMF.deafault}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get probanility of \eqn{x=1} conditioned by other variables.
}
\description{
probanility of \eqn{x=1} specified by \code{target}
}
\usage{
xprob_mNMF(...)

xprob_mNMF.default(varname,
                   vargroup,
                   V,
                   X,
                   Y,
                   W=NULL)

xprob_mNMF.formula(formula,
                   mNMFobj,
                   varname,
                   data = parent.frame(),
                   W=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{varname}{the target variable of X}
  \item{vargroup}{variable name of all the X}
  \item{mNMFobj}{output of \code{mNMF_vb} function}
  \item{X}{sparse design matrix}
  \item{V}{matrix of parameters}
    \item{formula}{
  an object of class "formula": a symbolic description of the model to be fitted.
}
  \item{data}{	
  an optional data frame, containing the variables in the formula.
}
\item{W}{Prior probability of \eqn{x=1}. If \code{NULL}, use mean of X}
}
\examples{
### see vignette
}
