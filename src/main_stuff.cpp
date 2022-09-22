#include "dont_compile/fastnono.cpp"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List mixed_2group0(int nnt, int nn, int n, int k1, int k2, int k,
                         Eigen::MatrixXd a, Eigen::VectorXd y, Eigen::VectorXd ss,
                         double sigy, double sig1) {
  fit_out fit = mixed_2group(nnt, nn, n, k1, k2, k, a, y, ss, sigy, sig1);
  return Rcpp::List::create(Rcpp::Named("means") = fit.means,
                            Rcpp::Named("sds") = fit.sds,
                            Rcpp::Named("cov") = fit.cov,
                            Rcpp::Named("time") = fit.time);
}
