mNMF_vb_sp <- function(y, X, N0, probX0, L, 
                       iter=1000,
                       a=0.5, b=0.01,
                       V=NULL,
                       display_progress=TRUE,
                       indices=NULL){
  stopifnot(grepl(".gCMatrix",class(X)))
  if(is.null(indices)){
    indices <- attr(X, "indices")
  }
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  out <- doVB_pois_sp2(N = length(y),
                       yv = y,
                       xi = X@i, xp = X@p,
                       varind = indices,
                       probX0 = probX0, N0 = N0,
                       D = X@Dim[2],
                       L = L,
                       iter = iter,
                       a=a, b=b,
                       V = V,
                       display_progress=display_progress)
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

mNMF_svb_sp <- function(y, X, N, probX0, L,
                        n_batches,
                        n_epochs,
                        lr_param,
                        lr_type,
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
  out <- doSVB_pois_sp2(N, y,
                        X@i, X@p, varind = indices,
                        probX0 = probX0, N0 = N - length(y),
                        D = X@Dim[2],
                        L=L, iter=n_epochs, a=a, b=b,
                        V=V,
                        bsize = n_batches,
                        lr_param=lr_param,
                        lr_type=lr_type,
                        display_progress=display_progress)
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

mNMF_svb <- function(y, X, L,
                     n_batches,
                     n_epochs,
                     lr_param,
                     lr_type,
                     a = 1, b=1,
                     V=NULL,
                     display_progress=TRUE,
                     indices=NULL){
  stopifnot(grepl(".gCMatrix",class(X)))
  if(is.null(indices)){
    indices <- attr(X, "indices")
  }
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  if(class(y) == "dsparseVector" | class(y) == "isparseVector"){
    y <- as(y, "sparseVector")
  }
  out <- doSVB_pois_sp(y@length, y@x, y@i-1L,
                       X@i, X@p, indices, X@Dim[2],
                       L=L, iter=n_epochs, a=a, b=b,
                       V=V,
                       bsize = n_batches,
                       lr_param=lr_param,
                       lr_type=lr_type,
                       display_progress=display_progress)
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