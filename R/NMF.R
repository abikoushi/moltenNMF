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
  Vini <- list(matrix(rgamma(dims[1]*rank, 1, 10), dims[1], rank),
                 t(apply(Y, 2, sample, size=rank)+0.1))
  }
  naY = is.na(Y)
  if(any(naY)){
    weight <- list( (Y@Dim[1] - rowSums(naY))/Y@Dim[1],
                    (Y@Dim[2] - colSums(naY))/Y@Dim[2] )
    nn_wch = which(!is.na(Y@x))
    out = doVB_pois_2D_ww(Vini,
                       y = Y@x[nn_wch],
                       rowi = Y@i[nn_wch],  coli = Y@j[nn_wch],
                       dims = dims,
                       L = rank, iter=iter, 
                       weight = weight,
                       a = prior_shape, b = prior_rate,
                       display_progress = display_progress)
  }else{
    out = doVB_pois_2D(Vini,
                       y = Y@x, rowi = Y@i,  coli = Y@j,
                       dims = dims,
                       L = rank, iter=iter, 
                       a = prior_shape, b = prior_rate,
                       display_progress = display_progress)
  
  }
  return(out)
}


NMF2D_svb <- function(Y, rank,
                      n_epochs, 
                      n_baches,
                      Vini = NULL,
                      lr_param = c(1, 0.8),
                      lr_type = "power",
                      dims=NULL,
                      weight = NULL,
                      subiter = 1,
                      prior_shape=1, prior_rate=1,
                      index_decrement = 0,
                      display_progress=TRUE,
                      useonlyone=FALSE){
  if(all(class(Y)!="dgTMatrix")){
    Y = as(Y, "TsparseMatrix")    
  }
  if(is.null(dims)){
    dims <- dim(Y)
  }
  if(is.null(Vini)){
    Vini <- list(matrix(rgamma(dims[1]*rank, 1, 10), dims[1], rank),
                 t(apply(Y, 2, sample, size=rank)+0.1))
  }
  n_baches <- min(n_baches, length(Y@x))
  if(useonlyone){
    # out = doVB_pois_s_2D_t1(
    #   Vinit,
    #   y = Y@x, rowi = Y@i,  coli = Y@j,
    #   L = rank,
    #   iter = n_epochs,
    #   subiter = subiter,
    #   a = prior_shape, b = prior_rate,
    #   N1 = length(Y@x),
    #   Nr = dims[1],
    #   Nc = dims[2],
    #   bsize = n_baches,
    #   lr_param = lr_param,
    #   lr_type = lr_type,
    #   display_progress = display_progress)
  }else{
    out = doVB_pois_s_2D(y = Y@x, rowi = Y@i,  coli = Y@j,
                         L = rank,
                         iter = n_epochs,
                         subiter = subiter,
                         a = prior_shape, b = prior_rate,
                         N1 = length(Y@x),
                         Nr = dims[1],
                         Nc = dims[2],
                         bsize = n_baches,
                         lr_param = lr_param,
                         lr_type = lr_type,
                         display_progress = display_progress)    
  }
  return(out)
}

###

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
