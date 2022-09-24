#' @export rpredictor_mNMF
rpredictor_mNMF <- function(X, np, alpha, beta){
  return(return(PoissonGamma_rng(X@Dim[1], np, X@i, X@p, alpha, 1/beta)))
}

#' @export rpredictor_mrNMF
rpredictor_mrNMF <- function(X, np, alpha, beta, tau){
  return(return(NegbinGamma_rng(X@Dim[1], np, X@i, X@p, alpha, 1/beta, tau)))
}
