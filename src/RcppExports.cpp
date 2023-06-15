// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pochhammer
NumericVector pochhammer(const NumericVector a_real, const NumericVector a_imag, const NumericVector q_real, const NumericVector q_imag, const NumericVector n);
RcppExport SEXP _queueR_pochhammer(SEXP a_realSEXP, SEXP a_imagSEXP, SEXP q_realSEXP, SEXP q_imagSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type a_real(a_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type a_imag(a_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_real(q_realSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type q_imag(q_imagSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(pochhammer(a_real, a_imag, q_real, q_imag, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_queueR_pochhammer", (DL_FUNC) &_queueR_pochhammer, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_queueR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
