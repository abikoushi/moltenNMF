\name{learning_rate}
\alias{learning_rate}
\title{Apply a Learning Rate Schedule}
\description{
Calculates the learning rate at a specified step using the specified schedule.
}
\usage{
learning_rate(x, lr_param, lr_type)
}
\arguments{
  \item{x}{Numeric vector of input values (e.g., iterations or time steps).}
  \item{lr_param}{Parameters used to define the learning rate schedule.}
  \item{lr_type}{Type of learning rate schedule ("exponential" or "constanat").}
}
\examples{
curve(learning_rate(x, lr_param = c(15,0.7), lr_type = "exponential"),
      xlim = c(1,10))

curve(learning_rate(x, lr_param = c(0.05), lr_type = "const"),
      xlim = c(1,10))
}