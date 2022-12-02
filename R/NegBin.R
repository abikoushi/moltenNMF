#' @export mrNMF_vb.default
mrNMF_vb.default <- function(y, X, L,
                             iter=1000, a=0.5, b=0.01,
                             indices=NULL){
  stopifnot(class(X)=="lgCMatrix")
  if(is.null(indices)){
    indices <- attr(X,"indices")
  }
  out <- doVB_negbin(y, X@i, X@p,indices , X@Dim[2],
                 L=L, iter=iter, a=a, b=b)
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

#' @export mrNMF_vb.formula
mrNMF_vb.formula <- function(formula,
                            data = parent.frame(),
                            L, iter=1000, a=0.5, b=0.01){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doVB_negbin(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                     L=L, iter=iter, a=a, b=b)
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  vname <- attr(terms.formula(formula), "term.labels")
  di <- diff(ind)
  if(length(vname) < length(di)){
    vname <- c("(Intercept)", vname)
  }
  out[["vargroup"]] <- vname[rep(1:length(di),di)] 
  return(out)
}

#' @export
mrNMF_vb <- function(...){
  UseMethod("mrNMF_vb")
}
