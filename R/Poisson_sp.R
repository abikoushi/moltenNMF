# mNMF_svb_batch <- function(y, X, L, iter,
#                            N = NULL,
#                            a = 1, b=1,
#                            V=NULL,
#                            display_progress=TRUE,
#                            indices=NULL){
#   stopifnot(grepl(".gCMatrix",class(X)))
#   if(is.null(indices)){
#     indices <- attr(X, "indices")
#   }
#   if(is.null(V)){
#     V <- rinitV(ncol(X), L)
#   }
#   if(!any(class(y) == "dsparseVector" | class(y) == "isparseVector")){
#     stopifnot(!is.null(N))
#   out <- doVB_pois_sp_skip(N, y,
#                            X@i, X@p, indices, X@Dim[2],
#                            L = L, iter = iter, a = a, b = b,
#                            V = V, 
#                            display_progress=display_progress)
#   }else{
#     out <- doVB_pois_sp_skip(y@length, y@x,
#                              X@i, X@p, indices, X@Dim[2],
#                              L = L, iter = iter, a = a, b = b,
#                              V = V, 
#                              display_progress=display_progress)
#   }
#   rownames(out$shape) <- colnames(X)
#   rownames(out$rate) <- colnames(X)
#   if(!is.null(attr(X, "term.labels"))){
#     vname <- attr(X,"term.labels")
#     di <- diff(attr(X,"indices"))
#     if(length(vname) < length(di)){
#       vname <- c("(Intercept)", vname)
#     }
#     out[["vargroup"]] <- vname[rep(1:length(di),di)]
#   }
#   return(out)
# }

mNMF_svb <- function(y, X, L,
                     n_batches,
                     n_epochs,
                     lr_param,
                     lr_type,
                     N = NULL,
                     a = 1, b=1,
                     V=NULL,
                     display_progress=TRUE,
                     indices=NULL){
  stopifnot(grepl(".gCMatrix",class(X)))
  if(is.null(indices)){
    indices <- attr(X, "indices")
  }
  if(is.null(V)){
    V <- rinitV(ncol(X), L)
  }
  if(!any(class(y) == "dsparseVector" | class(y) == "isparseVector")){
    stopifnot(!is.null(N))
    out <- doSVB_pois_sp_skip(N, y,
                              X@i, X@p,
                              varind = indices,
                              D = X@Dim[2],
                              L=L, iter=n_epochs, a=a, b=b,
                              V=V,
                              bsize = n_batches,
                              lr_param=lr_param,
                              lr_type=lr_type,
                              display_progress=display_progress)    
  }else{
    out <- doSVB_pois_sp_skip(y@length, y@x,
                         X@i, X@p, 
                         varind = indices,
                         D = X@Dim[2],
                         L=L, iter=n_epochs, a=a, b=b,
                         V=V,
                         bsize = n_batches,
                         lr_param=lr_param,
                         lr_type=lr_type,
                         display_progress=display_progress)
  }
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  if(!is.null(attr(X, "term.labels"))){
    vname <- attr(X,"term.labels")
    di <- diff(attr(X,"indices"))
    if(length(vname) < length(di)){
      vname <- c("(Intercept)", vname)
    }
    out[["vargroup"]] <- vname[rep(1:length(di),di)]
  }
  return(out)
}

