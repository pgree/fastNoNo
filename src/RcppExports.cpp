// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mixed_2group
Rcpp::List mixed_2group(int nnt, int nn, int n, int k1, int k2, int k, Eigen::MatrixXd a, Eigen::VectorXd y, Eigen::VectorXd ss, double sigy, double sig1);
RcppExport SEXP _fastNoNo_mixed_2group(SEXP nntSEXP, SEXP nnSEXP, SEXP nSEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP kSEXP, SEXP aSEXP, SEXP ySEXP, SEXP ssSEXP, SEXP sigySEXP, SEXP sig1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nnt(nntSEXP);
    Rcpp::traits::input_parameter< int >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< int >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type a(aSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< double >::type sigy(sigySEXP);
    Rcpp::traits::input_parameter< double >::type sig1(sig1SEXP);
    rcpp_result_gen = Rcpp::wrap(mixed_2group(nnt, nn, n, k1, k2, k, a, y, ss, sigy, sig1));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastNoNo_mixed_2group", (DL_FUNC) &_fastNoNo_mixed_2group, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastNoNo(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
