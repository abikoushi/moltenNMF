#' @export mNMF_vb.default
mNMF_vb.default <- function(y, X, L, pos_missing = NULL, iter=1000,  a=1,  b=1, eta=1){
  stopifnot(class(X)=="lgCMatrix")
  if(!is.null(pos_missing)){
    sumx <- Matrix::colSums(X)
    sum1mx <- Matrix::colSums(1-X)
    out <- doVB_pois_missing(y, X@i, X@p,
                             miss_row = pos_missing[,1],
                             miss_col = pos_missing[,2],
                             sumx = sumx,
                             #sum1mx = sum1mx,
                             varind = attr(X,"indices"),
                             dropK = attr(X,"dropK"),
                             D = X@Dim[2],
                             L = L, iter=iter, a=a, b=b, eta=eta)
  }else{
    out <- doVB_pois(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                     L=L, iter=iter, a=a, b=b) 
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

#' @export mNMF_vb.formula
mNMF_vb.formula <- function(formula,
                            data = parent.frame(),
                            L, iter=1000, a=1, b=1){
  dat <- model.frame(formula, data=data)
  X <- sparse_model_matrix_b(formula, data=data)
  y <- model.response(dat)
  M <- is.na(dat)
  if(any(M)){
    posM <-  which(M, arr.ind = TRUE)-1L
    out <- doVB_pois_missing(y, X@i, X@p,
                             miss_row = posM[,1],
                             miss_col = posM[,2],
                             sumx = colSums(X),
                             #sum1mx = colSums(1-X),
                             varind = attr(X,"indices"),
                             dropK = attr(X,"dropK"),
                             D = X@Dim[2],
                             L=L, iter=iter, a=a, b=b)
  }else{
    out <- doVB_pois(y, X@i, X@p, attr(X,"indices"), X@Dim[2],
                     L=L, iter=iter, a=a, b=b) 
  }
  rownames(out$shape) <- colnames(X)
  rownames(out$rate) <- colnames(X)
  out[["X"]] <- X
  ind <- attr(X,"indices")
  #fchar <- as.character(formula)
  #vname <- unlist(strsplit(fchar[length(fchar)]," [+] | [-] "))
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
