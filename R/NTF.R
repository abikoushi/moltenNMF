# multidimensional array
NTF_vb <- function(Y, X, rank, 
                   weight = NULL,
                   dims = NULL, 
                   iter = 100, dict = NULL,
                   prior_shape=1, prior_rate=1,
                   index_decrement = 1L,
                   display_progress=TRUE){
  if(is.null(dims)){
    dims <- apply(X, 2, function(x){ diff(range(x)) }) + 1L
  }
  if(is.null(weight)){
    if(ncol(X)==2L){
      out <- doVB_pois_2D(y = Y, 
                          rowi = X[,1] - index_decrement, 
                          coli = X[,2] - index_decrement,
                          dims = dims,
                          L = rank, iter=iter,
                          a = prior_shape, b = prior_rate,
                          display_progress = display_progress)
    }else{
      out <- doVB_pois_arr(y = Y, X = X - index_decrement, dims = dims,
                       L = rank, iter=iter,
                       a = prior_shape, b = prior_rate, 
                       display_progress = display_progress)
    }
  }else{
    out <- doVB_pois_w_arr(y = Y, X = X - index_decrement, dims = dims,
                           L = rank, iter=iter,
                           a = prior_shape, b = prior_rate,
                           display_progress = display_progress,
                           weight = weight)
  }
  
  names(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  if(!is.null(dict)){
    for(i in 1:length(out$shape)){
      row.names(out$shape[[i]]) <- dict[[i]]
    }
  }
  return(out)
}

NTF_svb <- function(Y, X,
                    rank,
                    n_epochs,
                    n_baches,
                    lr_param = c(1, 0.8),
                    lr_type = "exponential",
                    dims=NULL,
                    weight = NULL,
                    subiter = 1,
                    prior_shape=1, prior_rate=1,
                    index_decrement = 1L,
                    display_progress=TRUE){
  if(is.null(dims)){
    dims <- apply(X, 2, function(x)diff(range(x)))+1L
  }
  n_baches <- min(n_baches, length(Y))
  if(ncol(X)==2L){
    out = doVB_pois_s_2D(y = Y, rowi = X@i,  coli = X@j,
                         L = rank,
                         iter = n_epochs,
                         subiter = subiter,
                         a = prior_shape, b = prior_rate,
                         N1 = length(Y),
                         Nr = dims[1],
                         Nc = dims[2],
                         bsize = n_baches,
                         lr_param = lr_param,
                         lr_type = lr_type,
                         display_progress = display_progress)
  }else{
    out <- doVB_pois_s_arr(y = Y,
                           X = X - index_decrement,
                           dims = dims,
                           L = rank,
                           iter = n_epochs,
                           subiter = subiter,
                           a = prior_shape, b = prior_rate,
                           N1 = length(Y),
                           bsize = n_baches,
                           lr_param = lr_param,
                           lr_type = lr_type,
                           display_progress)
  }
  return(out)
}

