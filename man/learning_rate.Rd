\name{learning_rate}
\alias{learning_rate}
\title{Apply a Learning Rate Schedule}
\description{
  Applies a learning rate transformation to each element of a numeric vector using a specified schedule.
}
\usage{
learning_rate(x, lr_param, lr_type)
}
\arguments{
  \item{x}{Numeric vector of input values (e.g., iterations or time steps).}
  \item{lr_param}{Parameters used to define the learning rate schedule.}
  \item{lr_type}{Type of learning rate schedule ("exponential" or "constanat").}
}
\details{
  This function applies the helper function \code{check_lr} to each element of \code{x}, returning a vector of learning rates. The definition of \code{check_lr} is assumed to be available in the environment.
}
\value{
  A numeric vector of learning rates corresponding to input \code{x}.
}