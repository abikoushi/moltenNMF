#' @export rpredictor_mNMF
rpredictor_mNMF <- function(X, np, shape, rate){
  return(return(PoissonGamma_rng(X@Dim[1], np, X@i, X@p, shape, 1/rate)))
}

#' @export rpredictor_mrNMF
rpredictor_mrNMF <- function(X, np, alpha, shape, rate, precision){
  return(return(NegbinGamma_rng(X@Dim[1], np, X@i, X@p, shape, 1/rate, precision)))
}
