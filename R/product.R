#' @export product_m.formula
product_m.formula <- function(formula, data = parent.frame(), V, aggregate = TRUE){
  X <- sparse_onehot(formula, data=data)
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


product_array <- function(V, X){
  V0 <- V[[1]][X[,1],]
  if(ncol(X)>=2){
    for(i in 2:ncol(X)){
      V0 <- V0 * V[[i]][X[,i],]  
    }
  }
  rowSums(V0)
}
