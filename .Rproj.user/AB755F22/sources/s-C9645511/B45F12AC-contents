#' @export
sparse_model_matrix_b <- function(object,
                                  data = environment(object), ...){
  X <- Matrix::sparse.model.matrix(object, data = data,...)
  m <- attr(X,"assign")
  ind <- which(diff(m)!=0)
  Y <- methods::as(X, "lgCMatrix") # equivalent to X > 0
  attr(Y,"indices") <- c(0,ind,length(m))
  return(Y)
}

###
#' @export
prod_mPPCA.formula <- function(formula, data = parent.frame(), mu){
  X <- sparse_model_matrix_b(formula, data=data)
  return(myprod(X@Dim[1], X@i, X@p, mu))
}
#' @export
prod_mPPCA.dafault <- function(X, mu){
  return(myprod(X@Dim[1], X@i, X@p, mu))
}
#' @export
prod_mPPCA <- function(...){
  UseMethod("prod_mPPCA")
}

####
#' @export
mPPCA_Gibbs.dafault <- function(y, X, L, iter=2000, lambda=1, tau=1){
  out <- doGibbs(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}
#' @export
mPPCA_Gibbs.formula <- function(formula,
                                data = parent.frame(),
                                L, iter=2000, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doGibbs(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}
#' @export
mPPCA_Gibbs <- function(...){
  UseMethod("mPPCA_Gibbs")
}

###
#' @export
mPPCA_vb.dafault <- function(y, X, L, iter=200, lambda=1, tau=1){
  out <- doVB(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mean) <- colnames(X)
  rownames(out$sd) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}
#' @export
mPPCA_vb.formula <- function(formula,
                                data = parent.frame(),
                                L, iter=200, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doVB(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mean) <- colnames(X)
  rownames(out$sd) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}
#' @export
mPPCA_vb <- function(...){
  UseMethod("mPPCA_vb")
}
