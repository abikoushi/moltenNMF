meanV_array <- function(out, logarithm = FALSE){
  if(logarithm){
    V <- lapply(1:length(out$shape), function(i)sweep(digamma(out$shape[[i]]), 2, log(out$rate[i,]), "-"))
  }else{
    V <- lapply(1:length(out$shape), function(i)sweep(out$shape[[i]], 2, out$rate[i,], "/"))
  }
  names(V) <- names(out$shape)
  return(V)
}

meanV <- function(out, logarithm = FALSE){
  if(is.list(out$shape)){
    V = meanV_array(out = out, logarithm = logarithm)
  }
  if(is.matrix(out$shape)){
    if(logarithm){
      V = digamma(out$shape) - log(rate)
    }else{
      V = out$shape/out$rate      
    }
  }
  return(V)
}

rearrange_cols <- function(Vm, axis = 1L,
                           FUN = var, 
                           normalize = FALSE,
                           decreasing = TRUE){
  if(is.matrix(Vm)){
    if(normalize){
      Vm <- sweep(Vm, 1, rowSums(Vm), FUN = "/") 
    }
    ord = order(apply(Vm, 2, FUN = FUN), decreasing = decreasing)
    Vm = Vm[, ord]
  }
  if(is.list(Vm)){
    if(normalize){
      for(i in 1:length(Vm)){
        Vm[[i]] <- sweep(Vm[[i]], 1, rowSums(Vm[[i]]), FUN = "/")      
      }
    }
    ord = order(apply(Vm[[axis]], 2, FUN = FUN), decreasing = decreasing)
    for(i in 1:length(Vm)){
      Vm[[i]] = Vm[[i]][, ord]
    }    
  }
  return(Vm)
}
