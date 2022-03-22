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
