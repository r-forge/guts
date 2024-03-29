// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// guts_engine
void guts_engine(Rcpp::List gobj, Rcpp::NumericVector par, Rcpp::Nullable<Rcpp::NumericVector > z_dist);
RcppExport SEXP _GUTS_guts_engine(SEXP gobjSEXP, SEXP parSEXP, SEXP z_distSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type gobj(gobjSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector > >::type z_dist(z_distSEXP);
    guts_engine(gobj, par, z_dist);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GUTS_guts_engine", (DL_FUNC) &_GUTS_guts_engine, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_GUTS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
