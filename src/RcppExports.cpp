// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dcca
List dcca(arma::vec x, arma::vec y, int order, arma::ivec scales);
RcppExport SEXP _fractalRegression_dcca(SEXP xSEXP, SEXP ySEXP, SEXP orderSEXP, SEXP scalesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type scales(scalesSEXP);
    rcpp_result_gen = Rcpp::wrap(dcca(x, y, order, scales));
    return rcpp_result_gen;
END_RCPP
}
// dfa
List dfa(arma::vec x, int order, arma::uword verbose, arma::uvec scales, double scale_ratio);
RcppExport SEXP _fractalRegression_dfa(SEXP xSEXP, SEXP orderSEXP, SEXP verboseSEXP, SEXP scalesSEXP, SEXP scale_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< double >::type scale_ratio(scale_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(dfa(x, order, verbose, scales, scale_ratio));
    return rcpp_result_gen;
END_RCPP
}
// dlcca
arma::mat dlcca(arma::vec x, arma::vec y, int order, IntegerVector scales, int lags, char direction);
RcppExport SEXP _fractalRegression_dlcca(SEXP xSEXP, SEXP ySEXP, SEXP orderSEXP, SEXP scalesSEXP, SEXP lagsSEXP, SEXP directionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< char >::type direction(directionSEXP);
    rcpp_result_gen = Rcpp::wrap(dlcca(x, y, order, scales, lags, direction));
    return rcpp_result_gen;
END_RCPP
}
// poly_residuals
arma::vec poly_residuals(arma::vec yr, int m);
RcppExport SEXP _fractalRegression_poly_residuals(SEXP yrSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type yr(yrSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_residuals(yr, m));
    return rcpp_result_gen;
END_RCPP
}
// lm_c
arma::vec lm_c(arma::vec xs, arma::vec yr);
RcppExport SEXP _fractalRegression_lm_c(SEXP xsSEXP, SEXP yrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yr(yrSEXP);
    rcpp_result_gen = Rcpp::wrap(lm_c(xs, yr));
    return rcpp_result_gen;
END_RCPP
}
// seq_int
arma::uvec seq_int(arma::uword length);
RcppExport SEXP _fractalRegression_seq_int(SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_int(length));
    return rcpp_result_gen;
END_RCPP
}
// seq_int_range
arma::uvec seq_int_range(arma::uword start, arma::uword stop);
RcppExport SEXP _fractalRegression_seq_int_range(SEXP startSEXP, SEXP stopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword >::type start(startSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type stop(stopSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_int_range(start, stop));
    return rcpp_result_gen;
END_RCPP
}
// detrend_cov
arma::vec detrend_cov(arma::vec x, arma::vec y, int m);
RcppExport SEXP _fractalRegression_detrend_cov(SEXP xSEXP, SEXP ySEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(detrend_cov(x, y, m));
    return rcpp_result_gen;
END_RCPP
}
// mfdfa
List mfdfa(arma::vec x, arma::vec q, int order, arma::uvec scales, double scale_ratio);
RcppExport SEXP _fractalRegression_mfdfa(SEXP xSEXP, SEXP qSEXP, SEXP orderSEXP, SEXP scalesSEXP, SEXP scale_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< double >::type scale_ratio(scale_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(mfdfa(x, q, order, scales, scale_ratio));
    return rcpp_result_gen;
END_RCPP
}
// mfdfa_cj
List mfdfa_cj(arma::vec Timeseries, arma::vec qValues, arma::uvec scales);
RcppExport SEXP _fractalRegression_mfdfa_cj(SEXP TimeseriesSEXP, SEXP qValuesSEXP, SEXP scalesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Timeseries(TimeseriesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type qValues(qValuesSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type scales(scalesSEXP);
    rcpp_result_gen = Rcpp::wrap(mfdfa_cj(Timeseries, qValues, scales));
    return rcpp_result_gen;
END_RCPP
}
// win_sums
arma::vec win_sums(arma::vec x, uint window);
RcppExport SEXP _fractalRegression_win_sums(SEXP xSEXP, SEXP windowSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< uint >::type window(windowSEXP);
    rcpp_result_gen = Rcpp::wrap(win_sums(x, window));
    return rcpp_result_gen;
END_RCPP
}
// fitting
arma::vec fitting(arma::vec x, arma::vec y);
RcppExport SEXP _fractalRegression_fitting(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(fitting(x, y));
    return rcpp_result_gen;
END_RCPP
}
// mlra
List mlra(arma::vec x, arma::vec y, int order, IntegerVector scales, int lags, char direction);
RcppExport SEXP _fractalRegression_mlra(SEXP xSEXP, SEXP ySEXP, SEXP orderSEXP, SEXP scalesSEXP, SEXP lagsSEXP, SEXP directionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type scales(scalesSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< char >::type direction(directionSEXP);
    rcpp_result_gen = Rcpp::wrap(mlra(x, y, order, scales, lags, direction));
    return rcpp_result_gen;
END_RCPP
}
// lagn
arma::mat lagn(arma::vec x, int lags);
RcppExport SEXP _fractalRegression_lagn(SEXP xSEXP, SEXP lagsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type lags(lagsSEXP);
    rcpp_result_gen = Rcpp::wrap(lagn(x, lags));
    return rcpp_result_gen;
END_RCPP
}
// polyfit
arma::vec polyfit(arma::vec x, int order);
RcppExport SEXP _fractalRegression_polyfit(SEXP xSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(polyfit(x, order));
    return rcpp_result_gen;
END_RCPP
}
// detrend_var
arma::mat detrend_var(arma::mat X, int order);
RcppExport SEXP _fractalRegression_detrend_var(SEXP XSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(detrend_var(X, order));
    return rcpp_result_gen;
END_RCPP
}
// mra
List mra(arma::vec x, arma::vec y, int order, arma::ivec scales);
RcppExport SEXP _fractalRegression_mra(SEXP xSEXP, SEXP ySEXP, SEXP orderSEXP, SEXP scalesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< arma::ivec >::type scales(scalesSEXP);
    rcpp_result_gen = Rcpp::wrap(mra(x, y, order, scales));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fractalRegression_dcca", (DL_FUNC) &_fractalRegression_dcca, 4},
    {"_fractalRegression_dfa", (DL_FUNC) &_fractalRegression_dfa, 5},
    {"_fractalRegression_dlcca", (DL_FUNC) &_fractalRegression_dlcca, 6},
    {"_fractalRegression_poly_residuals", (DL_FUNC) &_fractalRegression_poly_residuals, 2},
    {"_fractalRegression_lm_c", (DL_FUNC) &_fractalRegression_lm_c, 2},
    {"_fractalRegression_seq_int", (DL_FUNC) &_fractalRegression_seq_int, 1},
    {"_fractalRegression_seq_int_range", (DL_FUNC) &_fractalRegression_seq_int_range, 2},
    {"_fractalRegression_detrend_cov", (DL_FUNC) &_fractalRegression_detrend_cov, 3},
    {"_fractalRegression_mfdfa", (DL_FUNC) &_fractalRegression_mfdfa, 5},
    {"_fractalRegression_mfdfa_cj", (DL_FUNC) &_fractalRegression_mfdfa_cj, 3},
    {"_fractalRegression_win_sums", (DL_FUNC) &_fractalRegression_win_sums, 2},
    {"_fractalRegression_fitting", (DL_FUNC) &_fractalRegression_fitting, 2},
    {"_fractalRegression_mlra", (DL_FUNC) &_fractalRegression_mlra, 6},
    {"_fractalRegression_lagn", (DL_FUNC) &_fractalRegression_lagn, 2},
    {"_fractalRegression_polyfit", (DL_FUNC) &_fractalRegression_polyfit, 2},
    {"_fractalRegression_detrend_var", (DL_FUNC) &_fractalRegression_detrend_var, 2},
    {"_fractalRegression_mra", (DL_FUNC) &_fractalRegression_mra, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_fractalRegression(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
