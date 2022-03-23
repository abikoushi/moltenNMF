#' @export mBMF_Gibbs.default
mBMF_Gibbs.default <- function(y, X, L, iter=2000, lambda=1, tau=1){
  out <- doGibbs_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                        L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}

#' @export mBMF_Gibbs.formula
mBMF_Gibbs.formula <- function(formula,
                               data = parent.frame(),
                               L, iter=2000, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doGibbs_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}

#' @export 
mBMF_Gibbs <- function(...){
  UseMethod("mBMF_Gibbs")
}
