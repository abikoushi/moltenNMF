#' @export prod_mNMF.formula
prod_mNMF.formula <- function(formula, data = parent.frame(), mu){
  X <- sparse_model_matrix_b(formula, data=data)
  return(myprod(X@Dim[1], X@i, X@p, mu))
}

#' @export prod_mNMF.default
prod_mNMF.default <- function(X, mu){
  return(myprod(X@Dim[1], X@i, X@p, mu))
}

#' @export
prod_mNMF <- function(...){
  UseMethod("prod_mNMF")
}
