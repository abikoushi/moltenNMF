#' @export prod_mPPCA.formula
prod_mPPCA.formula <- function(formula, data = parent.frame(), mu){
  X <- sparse_model_matrix_b(formula, data=data)
  return(myprod(X@Dim[1], X@i, X@p, mu))
}

#' @export prod_mPPCA.default
prod_mPPCA.default <- function(X, mu){
  return(myprod(X@Dim[1], X@i, X@p, mu))
}

#' @export
prod_mPPCA <- function(...){
  UseMethod("prod_mPPCA")
}
