// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// scr_nll_acoustic
double scr_nll_acoustic(NumericVector pars, NumericMatrix caps, NumericMatrix traps, NumericMatrix mask, NumericMatrix maskDists, NumericVector nCalls);
RcppExport SEXP _scr_scr_nll_acoustic(SEXP parsSEXP, SEXP capsSEXP, SEXP trapsSEXP, SEXP maskSEXP, SEXP maskDistsSEXP, SEXP nCallsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type caps(capsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type maskDists(maskDistsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nCalls(nCallsSEXP);
    rcpp_result_gen = Rcpp::wrap(scr_nll_acoustic(pars, caps, traps, mask, maskDists, nCalls));
    return rcpp_result_gen;
END_RCPP
}
// eucdist_nll
NumericMatrix eucdist_nll(NumericMatrix points, NumericMatrix traplocations);
RcppExport SEXP _scr_eucdist_nll(SEXP pointsSEXP, SEXP traplocationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type points(pointsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type traplocations(traplocationsSEXP);
    rcpp_result_gen = Rcpp::wrap(eucdist_nll(points, traplocations));
    return rcpp_result_gen;
END_RCPP
}
// scr_nll
double scr_nll(NumericVector pars, NumericMatrix caps, NumericMatrix traps, NumericMatrix mask, NumericMatrix maskDists);
RcppExport SEXP _scr_scr_nll(SEXP parsSEXP, SEXP capsSEXP, SEXP trapsSEXP, SEXP maskSEXP, SEXP maskDistsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type caps(capsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type traps(trapsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mask(maskSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type maskDists(maskDistsSEXP);
    rcpp_result_gen = Rcpp::wrap(scr_nll(pars, caps, traps, mask, maskDists));
    return rcpp_result_gen;
END_RCPP
}
// pointgen
NumericMatrix pointgen(int n, NumericVector xlim, NumericVector ylim);
RcppExport SEXP _scr_pointgen(SEXP nSEXP, SEXP xlimSEXP, SEXP ylimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xlim(xlimSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ylim(ylimSEXP);
    rcpp_result_gen = Rcpp::wrap(pointgen(n, xlim, ylim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scr_scr_nll_acoustic", (DL_FUNC) &_scr_scr_nll_acoustic, 6},
    {"_scr_eucdist_nll", (DL_FUNC) &_scr_eucdist_nll, 2},
    {"_scr_scr_nll", (DL_FUNC) &_scr_scr_nll, 5},
    {"_scr_pointgen", (DL_FUNC) &_scr_pointgen, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
