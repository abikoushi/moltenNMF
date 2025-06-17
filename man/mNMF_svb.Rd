\name{mNMF_svb}
\alias{mNMF_svb}
\alias{mNMF_svb_batch}
\title{
Estimate the model of molten NMF using stochastic variational Bayes
}
\description{
Performs mini-batch Stochastic Variational Bayes (SVB) inference for molten non-negative matrix factorization (NMF) under a Poisson likelihood on sparse input data.
}
\usage{
mNMF_svb(y, X, L,
         n_batches,
         n_epochs,
         lr_param,
         lr_type,
         N = NULL,
         a = 1, b = 1,
         V = NULL,
         display_progress = TRUE,
         indices = NULL)
         
mNMF_svb_batch(y, X, L, N, iter,
               a = 1, b=1,
               V=NULL,
               display_progress=TRUE,
               indices=NULL)
}
\arguments{
  \item{y}{A sparse vector of observed counts, either of class \code{dsparseVector} or \code{isparseVector}, or a dense numeric vector.}
  \item{X}{A sparse design matrix of class \code{dgCMatrix} or similar (\code{.gCMatrix}).}
  \item{L}{An integer specifying the number of latent components.}
  \item{n_batches}{Number of mini-batches per epoch.}
  \item{n_epochs}{Total number of epochs (iterations over the entire dataset).}
  \item{lr_param}{Learning rate parameters used in the SVB updates. Typically a scalar or a list.}
  \item{lr_type}{Type of learning rate schedule (\code{"poly"} or \code{"exp"}).}
  \item{N}{Total number of possible observations. Required if \code{y} is dense.}
  \item{a}{Shape parameter of the Gamma prior. Default is \code{1}.}
  \item{b}{Rate parameter of the Gamma prior. Default is \code{1}.}
  \item{V}{Optional initialization of the variational parameter matrix \code{V}. If \code{NULL}, random initialization is used.}
  \item{display_progress}{Logical; whether to display iteration progress. Default is \code{TRUE}.}
  \item{indices}{Optional grouping indices for variable blocks (used for structured sparsity). If \code{NULL}, extracted from \code{X}.}
}
\value{
A list containing:
  \item{shape}{Matrix of variational shape parameters for the factors.}
  \item{rate}{Matrix of variational rate parameters for the factors.}
  \item{vargroup}{(Optional) A character vector indicating the variable group each column belongs to, if term labels are available.}
}
\details{
This function is designed to efficiently perform variational inference for non-negative matrix factorization under a Poisson observation model using mini-batch stochastic optimization. It supports sparse matrix input and sparse response vectors, with grouped variable structures tracked via attributes of \code{X}.
}
\examples{
library(Matrix)

# Simulated sparse data
X <- rsparsematrix(1000, 20, density = 0.05)
y <- sample(0:3, 1000, replace = TRUE)
y_sparse <- as(y, "dsparseVector")

# Run stochastic variational inference
result <- mNMF_svb(y_sparse, X, L = 5,
                   n_batches = 10,
                   n_epochs = 50,
                   lr_param = 0.1,
                   lr_type = "poly",
                   a = 1, b = 1)
}
\seealso{
\code{\link{mNMF_vb}}
}
