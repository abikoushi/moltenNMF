logsumexp <- function(X){
  mx <- apply(X, 1, max)
  return(mx + log(Matrix::rowSums(exp(sweep(X,1,mx)))))
}

xprob <- function(V,target,W,X,Y){
  lvt <- log(V[target,,drop=FALSE])
  tlen <- nrow(lvt)
  prob <- matrix(0,length(Y),tlen)
  wt <- W[target]
  t0 <- X[,!target,drop=FALSE]%*%log(V[!target,,drop=FALSE])
  for(i in 1:tlen){
    t1 <- sweep(t0, 2, lvt[i,], "+")
    T1 <- Y*logsumexp(t1)-Matrix::rowSums(exp(t1))+log(wt[i])
    T0 <- Y*logsumexp(t0)-Matrix::rowSums(exp(t0))+log1p(-wt[i])
    prob[,i] <- 1/(1+exp(T0-T1))
  }
  return(prob)
}

#' @export xprob_mNMF.formula
xprob_mNMF.formula <- function(formula,
                               mNMFobj,
                               varname,
                               data = parent.frame(),
                               W=NULL){
  X <- sparse_model_matrix_b(formula, data)
  Y <- model.response(model.frame(formula, data=data))
  if(is.null(W)){
    W <- colMeans(X, na.rm = TRUE) 
  }
  V <- mNMFobj$shape/mNMFobj$rate
  target <- mNMFobj$vargroup==as.character(varname)
  return(xprob(V,target,W,X,Y))
}

#' @export xprob_mNMF.default
xprob_mNMF.default <- function(varname,
                               vargroup,
                               V,
                               X,
                               Y,
                               W=NULL){
  if(is.null(W)){
    W <- Matrix::colMeans(X, na.rm = TRUE) 
  }
  target <- vargroup==as.character(varname)
  return(xprob(V,target,W,X,Y))
}

#' @export
xprob_mNMF <- function(...){
  UseMethod("xprob_mNMF")
}
