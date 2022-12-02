#' @export mNMF_vb.default
mNMF_vb.default <- function(y, X, L,
                            iter=1000, a=0.5, b=0.01,
                            display_progress=TRUE,
                            indices=NULL){
  stopifnot(class(X)=="lgCMatrix")
  if(is.null(indices)){
        indices <- attr(X,"indices")
  }
  out <- doVB_pois(y, X@i, X@p, indices, X@Dim[2],
                     L=L, iter=iter, a=a, b=b,
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

#' @export mNMF_vb.formula
mNMF_vb.formula <- function(formula,
                            data = parent.frame(),
                            L, iter=1000, a=0.5, b=0.01,
                            display_progress=TRUE){
  dat <- model.frame(formula, data=data)
  X <- sparse_onehot(formula, data=data)
  y <- model.response(dat)
  out <- doVB_pois(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                     L=L, iter=iter, a=a, b=b,
                     display_progress=display_progress) 
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

#' @export
mNMF_vb <- function(...){
  UseMethod("mNMF_vb")
}
