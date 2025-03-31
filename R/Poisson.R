rinitV <- function(D, L){
  matrix(rexp(D*L),D,L)
}

rinitV2 <- function(D, L, a, b){
  matrix(rgamma(D*L,a,b),D,L)
}

mNMF_vb.default <- function(y, X, L,
                            iter=1000,
                            a=0.5, b=0.01,
                            V=NULL,
                            offset = NULL,
                            display_progress=TRUE,
                            indices=NULL){
  stopifnot(grepl(".gCMatrix",class(X)))
  if(is.null(indices)){
    indices <- attr(X, "indices")
  }
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  if(is.null(offset)){
    if(class(y) == "dsparseVector" | class(y) == "isparseVector"){
      out <- doVB_pois_sp(y@length, y@x, y@i-1L,
                          X@i, X@p, indices, X@Dim[2],
                          L=L, iter=iter, a=a, b=b,
                          V=V,
                          display_progress=display_progress)
    }else{
      out <- doVB_pois(y, X@i, X@p, indices, X@Dim[2],
                       L=L, iter=iter, a=a, b=b,
                       V=V,
                       display_progress=display_progress) 
    }
  }else{
    if(class(y) == "dsparseVector" | class(y) == "isparseVector"){
      out <- doVB_pois_offset_sp(y@length, y@x, y@i-1L,
                                 X@i, X@p, indices, X@Dim[2],
                                 L=L, tau = offset,
                                 iter=iter, a=a, b=b,
                                 V=V,
                                 display_progress=display_progress)
    }else{
      out <- doVB_pois_offset(y,
                              X@i, X@p, indices, X@Dim[2],
                              L=L, tau = offset,
                              iter=iter, a=a, b=b,
                              V=V,
                              display_progress=display_progress)
    }
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

mNMF_vb.formula <- function(formula,
                            data = parent.frame(),
                            L, iter=1000, a=0.5, b=0.01,
                            V=NULL,
                            offset = NULL,
                            display_progress=TRUE){
  dat <- model.frame(formula, data=data)
  X <- sparse_onehot(formula, data=data)
  indices <- attr(X,"indices")
  y <- model.response(dat)
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  if(is.null(offset)){
      out <- doVB_pois(y, X@i, X@p, indices, X@Dim[2],
                       L=L, iter=iter, a=a, b=b,
                       V=V,
                       display_progress=display_progress) 
  }else{
      out <- doVB_pois_offset(y,
                              X@i, X@p, indices, X@Dim[2],
                              L=L, tau = offset,
                              iter=iter, a=a, b=b,
                              V=V,
                              display_progress=display_progress)
  }
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  out[["X"]] <- X
  ind <- attr(X,"indices")
  vname <- attr(terms.formula(formula), "term.labels")
  di <- diff(ind)
  if(length(vname) < length(di)){
    vname <- c("(Intercept)", vname)
  }
  out[["vargroup"]] <- vname[rep(1:length(di),di)] 
  return(out)
}

mNMF_vb <- function(...){
  UseMethod("mNMF_vb")
}


mrNMF_vb.default <- function(y, X, L,
                             iter=1000, a=0.5, b=0.01,
                             V=NULL,
                             indices=NULL,
                             display_progress=TRUE){
  stopifnot(grepl(".gCMatrix",class(X)))
  if(is.null(indices)){
    indices <- attr(X,"indices")
  }
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  out <- doVB_negbin(y, X@i, X@p,indices , X@Dim[2],
                     L=L, iter=iter, a=a, b=b,
                     V=V,
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

mrNMF_vb.formula <- function(formula,
                             data = parent.frame(),
                             L, iter=1000, a=0.5, b=0.01,
                             V=NULL){
  if(is.null(V)){
    V <- rinitV(ncol(X),L)
  }
  X <- sparse_onehot(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doVB_negbin(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                     L=L, iter=iter, a=a, b=b,
                     V=V,
                     display_progress=display_progress)
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  vname <- attr(terms.formula(formula), "term.labels")
  out[["X"]] <- X
  ind <- attr(X,"indices")
  di <- diff(ind)
  if(length(vname) < length(di)){
    vname <- c("(Intercept)", vname)
  }
  out[["vargroup"]] <- vname[rep(1:length(di),di)] 
  return(out)
}

mrNMF_vb <- function(...){
  UseMethod("mrNMF_vb")
}

####

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