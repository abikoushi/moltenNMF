# 2D array
NMF2D_vb <- function(Y, rank,
                     iter=100,
                     Vini = NULL,
                     dims=NULL, 
                     prior_shape = 1, prior_rate = 1,
                     display_progress=TRUE){
  if(all(class(Y)!="dgTMatrix")){
    Y = as(Y, "TsparseMatrix") 
  }
  if(is.null(dims)){
    dims <- dim(Y)
  }
  if(is.null(Vini)){
  #   Vini <- list(t(apply(Y, 1, sample, size=rank)+0.1),
  #               t(apply(Y, 2, sample, size=rank)+0.1))
    # Vini <- list(matrix(sample(Y@x, size=dims[1]*rank, replace = TRUE), dims[1], rank),
    #              matrix(sample(Y@x, size=dims[2]*rank, replace = TRUE), dims[2], rank))
    Vini <- list(matrix(rgamma(dims[1]*rank, 1, 10), dims[1], rank),
                 t(apply(Y, 2, sample, size=rank)+0.1))
  }
  doVB_pois_2D(Vini,
               y = Y@x, rowi = Y@i,  coli = Y@j,
               dims = dims,
               L = rank, iter=iter,
               a = prior_shape, b = prior_rate,
               display_progress = display_progress)
}

meanV_2D <- function(out){
  V <- lapply(1:length(out$shape), function(i)sweep(out$shape[[i]], 2, out$rate[i,],"/"))
  names(V) <- names(out$shape)
  return(V)
}

rearrange_cols <- function(Vm, axis=1, FUN = mean, normalize = TRUE, decreasing = TRUE){
  if(normalize){
    for(i in 1:length(Vm)){
      Vm[[i]] <- sweep(Vm[[i]], 1, rowSums(Vm[[i]]), FUN = "/")      
    }
  }
  ord = order(apply(Vm[[axis]], 2, FUN = FUN), decreasing = decreasing)
  for(i in 1:length(Vm)){
    Vm[[i]] <- Vm[[i]][, ord]
  }
  return(Vm)
}
