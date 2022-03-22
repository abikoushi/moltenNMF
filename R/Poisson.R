#' @export mNMF_vb.default
mNMF_vb.default <- function(y, X, L, iter=1000,  a=1,  b=1){
  out <- doVB_pois(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, a=a, b=b)
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}

#' @export mNMF_vb.formula
mNMF_vb.formula <- function(formula,
                                data = parent.frame(),
                                L, iter=1000, a=1, b=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doVB_pois(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, a=a, b=b)
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  out[["indices"]]=attr(X,"indices")
  return(out)
}

#' @export
mNMF_vb <- function(...){
  UseMethod("mNMF_vb")
}
