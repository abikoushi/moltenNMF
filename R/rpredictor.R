#' @export rpredictor_mNMF
rpredictor_mNMF <- function(X, np, shape, rate){
  PoissonGamma_rng(X@Dim[1], np, X@i, X@p, shape, rate)
}

#' @export rpredictor_mrNMF
rpredictor_mrNMF <- function(X, np, shape, rate, precision){
  NegbinGamma_rng(X@Dim[1], np, X@i, X@p, shape, rate, precision)
}

#' @export dpredictor_mNMF
dpredictor_mNMF <- function(y, X, np, shape, rate){
  Poisson_lp(y, np, X@i, X@p, shape, rate)
}

#' @export dpredictor_mNMF
dpredictor_mrNMF <- function(y, X, np, shape, rate, precision){
  NegBin_lp(y, np, X@i, X@p, shape, rate, precision)
}
