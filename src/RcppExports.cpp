// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// findmiss
IntegerVector findmiss(NumericVector x);
RcppExport SEXP _linelistBayes_findmiss(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(findmiss(x));
    return rcpp_result_gen;
END_RCPP
}
// get_mu_vec
NumericVector get_mu_vec(NumericMatrix x12, NumericVector beta);
RcppExport SEXP _linelistBayes_get_mu_vec(SEXP x12SEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x12(x12SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_mu_vec(x12, beta));
    return rcpp_result_gen;
END_RCPP
}
// dummy
NumericMatrix dummy(IntegerVector week, IntegerVector weekend);
RcppExport SEXP _linelistBayes_dummy(SEXP weekSEXP, SEXP weekendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type week(weekSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weekend(weekendSEXP);
    rcpp_result_gen = Rcpp::wrap(dummy(week, weekend));
    return rcpp_result_gen;
END_RCPP
}
// logLikNB
double logLikNB(NumericVector delay_vec, NumericMatrix x12, NumericVector disp, NumericVector betaplus, int maxdelay);
RcppExport SEXP _linelistBayes_logLikNB(SEXP delay_vecSEXP, SEXP x12SEXP, SEXP dispSEXP, SEXP betaplusSEXP, SEXP maxdelaySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type delay_vec(delay_vecSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x12(x12SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type disp(dispSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type betaplus(betaplusSEXP);
    Rcpp::traits::input_parameter< int >::type maxdelay(maxdelaySEXP);
    rcpp_result_gen = Rcpp::wrap(logLikNB(delay_vec, x12, disp, betaplus, maxdelay));
    return rcpp_result_gen;
END_RCPP
}
// lambda
NumericVector lambda(NumericVector curve, NumericVector si);
RcppExport SEXP _linelistBayes_lambda(SEXP curveSEXP, SEXP siSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type curve(curveSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type si(siSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda(curve, si));
    return rcpp_result_gen;
END_RCPP
}
// getr
NumericVector getr(NumericVector curve, NumericVector si, int size);
RcppExport SEXP _linelistBayes_getr(SEXP curveSEXP, SEXP siSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type curve(curveSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type si(siSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(getr(curve, si, size));
    return rcpp_result_gen;
END_RCPP
}
// backnow_cm
List backnow_cm(NumericVector outcome, NumericVector days, IntegerVector week, IntegerVector weekend, int iter, double sigma, int maxdelay, NumericVector si, int size, int workerID, int printProgress, Nullable<int> cd);
RcppExport SEXP _linelistBayes_backnow_cm(SEXP outcomeSEXP, SEXP daysSEXP, SEXP weekSEXP, SEXP weekendSEXP, SEXP iterSEXP, SEXP sigmaSEXP, SEXP maxdelaySEXP, SEXP siSEXP, SEXP sizeSEXP, SEXP workerIDSEXP, SEXP printProgressSEXP, SEXP cdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type outcome(outcomeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type days(daysSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type week(weekSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weekend(weekendSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type maxdelay(maxdelaySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type si(siSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type workerID(workerIDSEXP);
    Rcpp::traits::input_parameter< int >::type printProgress(printProgressSEXP);
    Rcpp::traits::input_parameter< Nullable<int> >::type cd(cdSEXP);
    rcpp_result_gen = Rcpp::wrap(backnow_cm(outcome, days, week, weekend, iter, sigma, maxdelay, si, size, workerID, printProgress, cd));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_linelistBayes_findmiss", (DL_FUNC) &_linelistBayes_findmiss, 1},
    {"_linelistBayes_get_mu_vec", (DL_FUNC) &_linelistBayes_get_mu_vec, 2},
    {"_linelistBayes_dummy", (DL_FUNC) &_linelistBayes_dummy, 2},
    {"_linelistBayes_logLikNB", (DL_FUNC) &_linelistBayes_logLikNB, 5},
    {"_linelistBayes_lambda", (DL_FUNC) &_linelistBayes_lambda, 2},
    {"_linelistBayes_getr", (DL_FUNC) &_linelistBayes_getr, 3},
    {"_linelistBayes_backnow_cm", (DL_FUNC) &_linelistBayes_backnow_cm, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_linelistBayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
