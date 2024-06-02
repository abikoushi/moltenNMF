// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// doVB_negbin
List doVB_negbin(const arma::vec& y, const arma::uvec& xi, const arma::uvec& xp, const arma::uvec& varind, const int& D, const int& L, const int& iter, const double& a, const double& b, arma::mat& V, const bool& display_progress);
RcppExport SEXP _moltenNMF_doVB_negbin(SEXP ySEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP, SEXP VSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< const int& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_negbin(y, xi, xp, varind, D, L, iter, a, b, V, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// doVB_pois
List doVB_pois(const arma::vec& y, const arma::uvec& xi, const arma::uvec& xp, const arma::uvec& varind, const int& D, const int& L, const int& iter, const double& a, const double& b, arma::mat& V, const bool& display_progress);
RcppExport SEXP _moltenNMF_doVB_pois(SEXP ySEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP, SEXP VSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< const int& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_pois(y, xi, xp, varind, D, L, iter, a, b, V, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// doVB_pois_sp
List doVB_pois_sp(const arma::vec& yv, const arma::uvec& yi, const int& N, const arma::uvec& xi, const arma::uvec& xp, const arma::uvec& varind, const int& D, const int& L, const int& iter, const double& a, const double& b, arma::mat& V, const bool& display_progress);
RcppExport SEXP _moltenNMF_doVB_pois_sp(SEXP yvSEXP, SEXP yiSEXP, SEXP NSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP, SEXP VSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type yv(yvSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type yi(yiSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< const int& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_pois_sp(yv, yi, N, xi, xp, varind, D, L, iter, a, b, V, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// doVB_pois_offset
List doVB_pois_offset(const arma::vec& y, const arma::uvec& xi, const arma::uvec& xp, const arma::uvec& varind, const int& D, const int& L, const arma::vec& tau, const int& iter, const double& a, const double& b, arma::mat& V, const bool& display_progress);
RcppExport SEXP _moltenNMF_doVB_pois_offset(SEXP ySEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP tauSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP, SEXP VSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< const int& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const int& >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const bool& >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_pois_offset(y, xi, xp, varind, D, L, tau, iter, a, b, V, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// myprod
arma::mat myprod(const int& N, const arma::uvec& xi, const arma::uvec& xp, const arma::mat& lam);
RcppExport SEXP _moltenNMF_myprod(SEXP NSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(myprod(N, xi, xp, lam));
    return rcpp_result_gen;
END_RCPP
}
// summyprod
arma::vec summyprod(const int& n, const arma::uvec& xi, const arma::uvec& xp, const arma::mat& lam);
RcppExport SEXP _moltenNMF_summyprod(SEXP nSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(summyprod(n, xi, xp, lam));
    return rcpp_result_gen;
END_RCPP
}
// myprod_r
arma::mat myprod_r(const int& N, const arma::uvec& xj, const arma::uvec& xp, const arma::mat& lam);
RcppExport SEXP _moltenNMF_myprod_r(SEXP NSEXP, SEXP xjSEXP, SEXP xpSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(myprod_r(N, xj, xp, lam));
    return rcpp_result_gen;
END_RCPP
}
// myprod_r_i
arma::mat myprod_r_i(const int& N, const arma::uvec& xj, const arma::uvec& xp, const arma::uvec id, const arma::mat& lam);
RcppExport SEXP _moltenNMF_myprod_r_i(SEXP NSEXP, SEXP xjSEXP, SEXP xpSEXP, SEXP idSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xj(xjSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type id(idSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(myprod_r_i(N, xj, xp, id, lam));
    return rcpp_result_gen;
END_RCPP
}
// NegBin_lp
arma::mat NegBin_lp(const arma::vec& y, const int& np, const arma::uvec& xi, const arma::uvec& xp, const arma::mat& alpha, const arma::mat& beta, const double& tau);
RcppExport SEXP _moltenNMF_NegBin_lp(SEXP ySEXP, SEXP npSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type np(npSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(NegBin_lp(y, np, xi, xp, alpha, beta, tau));
    return rcpp_result_gen;
END_RCPP
}
// Poisson_lp
arma::mat Poisson_lp(const arma::vec& y, const int& np, const arma::uvec& xi, const arma::uvec& xp, const arma::mat& alpha, const arma::mat& beta);
RcppExport SEXP _moltenNMF_Poisson_lp(SEXP ySEXP, SEXP npSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type np(npSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(Poisson_lp(y, np, xi, xp, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// PoissonGamma_rng
arma::mat PoissonGamma_rng(int N, int np, arma::uvec xi, arma::uvec xp, arma::mat alpha, arma::mat beta);
RcppExport SEXP _moltenNMF_PoissonGamma_rng(SEXP NSEXP, SEXP npSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type np(npSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(PoissonGamma_rng(N, np, xi, xp, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// NegbinGamma_rng
arma::mat NegbinGamma_rng(int N, int np, arma::uvec xi, arma::uvec xp, arma::mat alpha, arma::mat beta, double tau);
RcppExport SEXP _moltenNMF_NegbinGamma_rng(SEXP NSEXP, SEXP npSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type np(npSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(NegbinGamma_rng(N, np, xi, xp, alpha, beta, tau));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_moltenNMF_doVB_negbin", (DL_FUNC) &_moltenNMF_doVB_negbin, 11},
    {"_moltenNMF_doVB_pois", (DL_FUNC) &_moltenNMF_doVB_pois, 11},
    {"_moltenNMF_doVB_pois_sp", (DL_FUNC) &_moltenNMF_doVB_pois_sp, 13},
    {"_moltenNMF_doVB_pois_offset", (DL_FUNC) &_moltenNMF_doVB_pois_offset, 12},
    {"_moltenNMF_myprod", (DL_FUNC) &_moltenNMF_myprod, 4},
    {"_moltenNMF_summyprod", (DL_FUNC) &_moltenNMF_summyprod, 4},
    {"_moltenNMF_myprod_r", (DL_FUNC) &_moltenNMF_myprod_r, 4},
    {"_moltenNMF_myprod_r_i", (DL_FUNC) &_moltenNMF_myprod_r_i, 5},
    {"_moltenNMF_NegBin_lp", (DL_FUNC) &_moltenNMF_NegBin_lp, 7},
    {"_moltenNMF_Poisson_lp", (DL_FUNC) &_moltenNMF_Poisson_lp, 6},
    {"_moltenNMF_PoissonGamma_rng", (DL_FUNC) &_moltenNMF_PoissonGamma_rng, 6},
    {"_moltenNMF_NegbinGamma_rng", (DL_FUNC) &_moltenNMF_NegbinGamma_rng, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_moltenNMF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
