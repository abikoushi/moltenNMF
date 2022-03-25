#' @export mPPCA_Gibbs.default
mPPCA_Gibbs.default <- function(y, X, L, iter=2000, lambda=1, tau=1){
  stopifnot(class(X)=="lgCMatrix")
  out <- doGibbs(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  if(!is.null(attr(X,"indices"))){
    out[["indices"]]=attr(X,"indices") 
  }
  return(out)
}

#' @export mPPCA_Gibbs.formula
mPPCA_Gibbs.formula <- function(formula,
                                data = parent.frame(),
                                L, iter=2000, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doGibbs(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  ind <- attr(X,"indices")
  out[["indices"]] <- ind
  fchar <- as.character(formula)
  vname <- unlist(strsplit(fchar[length(fchar)]," [+] | [-] "))
  out[["vargroup"]] <- vname[rep(1:(length(ind)-1), diff(ind))]
  return(out)
}

#' @export 
mPPCA_Gibbs <- function(...){
  UseMethod("mPPCA_Gibbs")
}

#' @export mPPCA_vb.default
mPPCA_vb.default <- function(y, X, L, iter=200, lambda=1, tau=1){
  stopifnot(class(X)=="lgCMatrix")
  out <- doVB(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mean) <- colnames(X)
  rownames(out$sd) <- colnames(X)
  if(!is.null(attr(X,"indices"))){
    out[["indices"]]=attr(X,"indices") 
  }
  return(out)
}

#' @export mPPCA_vb.formula
mPPCA_vb.formula <- function(formula,
                                data = parent.frame(),
                                L, iter=200, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doVB(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                 L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mean) <- colnames(X)
  rownames(out$sd) <- colnames(X)
  ind <- attr(X,"indices")
  out[["indices"]] <- ind
  fchar <- as.character(formula)
  vname <- unlist(strsplit(fchar[length(fchar)]," [+] | [-] "))
  out[["vargroup"]] <- vname[rep(1:(length(ind)-1), diff(ind))]
  return(out)
}

#' @export
mPPCA_vb <- function(...){
  UseMethod("mPPCA_vb")
}
