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
List doVB_negbin(arma::vec y, arma::uvec xi, arma::uvec xp, arma::uvec varind, int D, int L, int iter, double a, double b);
RcppExport SEXP _moltenNMF_doVB_negbin(SEXP ySEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_negbin(y, xi, xp, varind, D, L, iter, a, b));
    return rcpp_result_gen;
END_RCPP
}
// doVB_pois
List doVB_pois(arma::vec y, arma::uvec xi, arma::uvec xp, arma::uvec varind, int D, int L, int iter, double a, double b);
RcppExport SEXP _moltenNMF_doVB_pois(SEXP ySEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP varindSEXP, SEXP DSEXP, SEXP LSEXP, SEXP iterSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type varind(varindSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(doVB_pois(y, xi, xp, varind, D, L, iter, a, b));
    return rcpp_result_gen;
END_RCPP
}
// myprod
arma::mat myprod(int n, arma::uvec xi, arma::uvec xp, arma::mat lam);
RcppExport SEXP _moltenNMF_myprod(SEXP nSEXP, SEXP xiSEXP, SEXP xpSEXP, SEXP lamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type lam(lamSEXP);
    rcpp_result_gen = Rcpp::wrap(myprod(n, xi, xp, lam));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_moltenNMF_doVB_negbin", (DL_FUNC) &_moltenNMF_doVB_negbin, 9},
    {"_moltenNMF_doVB_pois", (DL_FUNC) &_moltenNMF_doVB_pois, 9},
    {"_moltenNMF_myprod", (DL_FUNC) &_moltenNMF_myprod, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_moltenNMF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
