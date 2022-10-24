#' @export product_m.formula
product_m.formula <- function(formula, data = parent.frame(), V, aggregate = TRUE){
  X <- sparse_model_matrix_b(formula, data=data)
  if(aggregate){
    res <- as.vector(rowSums(myprod(nrow(X), X@i, X@p, V)))
  }else{
    res <- myprod(nrow(X), X@i, X@p, V)
  }
  return(res)
}

#' @export product_m.default
product_m.default <- function(X, V, aggregate = TRUE){
  if(aggregate){
    res <- as.vector(rowSums(myprod(nrow(X), X@i, X@p, V)))
  }else{
    res <- myprod(nrow(X), X@i, X@p, V)
  }
  return(res)
}

#' @export
product_m <- function(...){
  UseMethod("product_m")
}
