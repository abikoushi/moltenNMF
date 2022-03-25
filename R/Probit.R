#' @export mBMF_Gibbs.default
mBMF_Gibbs.default <- function(y, X, L, iter=2000, lambda=1, tau=1){
  stopifnot(class(X)=="lgCMatrix")
  out <- doGibbs_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                        L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mu) <- colnames(X)
  if(!is.null(attr(X,"indices"))){
    out[["indices"]]=attr(X,"indices") 
  }
  return(out)
}

#' @export mBMF_Gibbs.formula
mBMF_Gibbs.formula <- function(formula,
                               data = parent.frame(),
                               L, iter=2000, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doGibbs_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
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
mBMF_Gibbs <- function(...){
  UseMethod("mBMF_Gibbs")
}

###

#' @export mBMF_mcvb.default
mBMF_mcvb.default <- function(y, X, L, iter=1000, lambda=1, tau=1){
  stopifnot(class(X)=="lgCMatrix")
  out <- doMCVB_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                        L=L, iter=iter, lambda=lambda, tau=tau)
  rownames(out$mean) <- colnames(X)
  rownames(out$sd) <- colnames(X)
  if(!is.null(attr(X,"indices"))){
    out[["indices"]]=attr(X,"indices") 
  }
  return(out)
}

#' @export mBMF_mcvb.formula
mBMF_mcvb.formula <- function(formula,
                             data = parent.frame(),
                             L, iter=1000, lambda=1, tau=1){
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(model.frame(formula, data=data))
  out <- doMCVB_probit(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
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
mBMF_mcvb <- function(...){
  UseMethod("mBMF_mcvb")
}
