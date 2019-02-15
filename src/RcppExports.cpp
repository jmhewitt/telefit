// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#define EIGEN_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// r_maternCov
arma::mat r_maternCov(arma::mat dist, double scale, double range, double smoothness, double nugget);
RcppExport SEXP _telefit_r_maternCov(SEXP distSEXP, SEXP scaleSEXP, SEXP rangeSEXP, SEXP smoothnessSEXP, SEXP nuggetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< double >::type smoothness(smoothnessSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    rcpp_result_gen = Rcpp::wrap(r_maternCov(dist, scale, range, smoothness, nugget));
    return rcpp_result_gen;
END_RCPP
}
// r_maternArray
arma::vec r_maternArray(arma::vec dist, double scale, double range, double smoothness, double nugget);
RcppExport SEXP _telefit_r_maternArray(SEXP distSEXP, SEXP scaleSEXP, SEXP rangeSEXP, SEXP smoothnessSEXP, SEXP nuggetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< double >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< double >::type smoothness(smoothnessSEXP);
    Rcpp::traits::input_parameter< double >::type nugget(nuggetSEXP);
    rcpp_result_gen = Rcpp::wrap(r_maternArray(dist, scale, range, smoothness, nugget));
    return rcpp_result_gen;
END_RCPP
}
// r_mc2_rinvwishart
arma::mat r_mc2_rinvwishart(arma::mat V, double n);
RcppExport SEXP _telefit_r_mc2_rinvwishart(SEXP VSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(r_mc2_rinvwishart(V, n));
    return rcpp_result_gen;
END_RCPP
}
// r_mvrnorm_postKron
arma::mat r_mvrnorm_postKron(arma::vec y, arma::mat A, arma::mat B, int nSamples, bool precision);
RcppExport SEXP _telefit_r_mvrnorm_postKron(SEXP ySEXP, SEXP ASEXP, SEXP BSEXP, SEXP nSamplesSEXP, SEXP precisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type precision(precisionSEXP);
    rcpp_result_gen = Rcpp::wrap(r_mvrnorm_postKron(y, A, B, nSamples, precision));
    return rcpp_result_gen;
END_RCPP
}
// r_mvrnorm_post
arma::mat r_mvrnorm_post(arma::vec y, arma::mat Sigma, int nSamples, bool precision);
RcppExport SEXP _telefit_r_mvrnorm_post(SEXP ySEXP, SEXP SigmaSEXP, SEXP nSamplesSEXP, SEXP precisionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< int >::type nSamples(nSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type precision(precisionSEXP);
    rcpp_result_gen = Rcpp::wrap(r_mvrnorm_post(y, Sigma, nSamples, precision));
    return rcpp_result_gen;
END_RCPP
}
// r_qintnorm
arma::vec r_qintnorm(arma::vec breaks, double mu, double sigma);
RcppExport SEXP _telefit_r_qintnorm(SEXP breaksSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type breaks(breaksSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(r_qintnorm(breaks, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dtest
NumericVector dtest(NumericVector x, int m, int n, int k, Eigen::SparseMatrix< double > R, double q, double ldetR, Eigen::MatrixXd Sigma);
RcppExport SEXP _telefit_dtest(SEXP xSEXP, SEXP mSEXP, SEXP nSEXP, SEXP kSEXP, SEXP RSEXP, SEXP qSEXP, SEXP ldetRSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix< double > >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type ldetR(ldetRSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dtest(x, m, n, k, R, q, ldetR, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// test_gmrf_approx
List test_gmrf_approx(NumericVector y, NumericVector x0);
RcppExport SEXP _telefit_test_gmrf_approx(SEXP ySEXP, SEXP x0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type x0(x0SEXP);
    rcpp_result_gen = Rcpp::wrap(test_gmrf_approx(y, x0));
    return rcpp_result_gen;
END_RCPP
}
// test_ll
NumericVector test_ll(NumericVector y, NumericVector lambda);
RcppExport SEXP _telefit_test_ll(SEXP ySEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ll(y, lambda));
    return rcpp_result_gen;
END_RCPP
}
// r_dgemkmm
arma::mat r_dgemkmm(arma::mat A, arma::mat B, arma::mat C);
RcppExport SEXP _telefit_r_dgemkmm(SEXP ASEXP, SEXP BSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(r_dgemkmm(A, B, C));
    return rcpp_result_gen;
END_RCPP
}
// r_rwishart
arma::mat r_rwishart(arma::mat V, int n);
RcppExport SEXP _telefit_r_rwishart(SEXP VSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(r_rwishart(V, n));
    return rcpp_result_gen;
END_RCPP
}
// r_rinvwishart
arma::mat r_rinvwishart(arma::mat V, int n);
RcppExport SEXP _telefit_r_rinvwishart(SEXP VSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(r_rinvwishart(V, n));
    return rcpp_result_gen;
END_RCPP
}
// r_dgeikmm
arma::mat r_dgeikmm(int N, arma::mat A, arma::mat B);
RcppExport SEXP _telefit_r_dgeikmm(SEXP NSEXP, SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(r_dgeikmm(N, A, B));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP r_ll(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP r_stpcomposition(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP r_stpfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP r_svcfit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP r_svcpredict(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_telefit_r_maternCov", (DL_FUNC) &_telefit_r_maternCov, 5},
    {"_telefit_r_maternArray", (DL_FUNC) &_telefit_r_maternArray, 5},
    {"_telefit_r_mc2_rinvwishart", (DL_FUNC) &_telefit_r_mc2_rinvwishart, 2},
    {"_telefit_r_mvrnorm_postKron", (DL_FUNC) &_telefit_r_mvrnorm_postKron, 5},
    {"_telefit_r_mvrnorm_post", (DL_FUNC) &_telefit_r_mvrnorm_post, 4},
    {"_telefit_r_qintnorm", (DL_FUNC) &_telefit_r_qintnorm, 3},
    {"_telefit_dtest", (DL_FUNC) &_telefit_dtest, 8},
    {"_telefit_test_gmrf_approx", (DL_FUNC) &_telefit_test_gmrf_approx, 2},
    {"_telefit_test_ll", (DL_FUNC) &_telefit_test_ll, 2},
    {"_telefit_r_dgemkmm", (DL_FUNC) &_telefit_r_dgemkmm, 3},
    {"_telefit_r_rwishart", (DL_FUNC) &_telefit_r_rwishart, 2},
    {"_telefit_r_rinvwishart", (DL_FUNC) &_telefit_r_rinvwishart, 2},
    {"_telefit_r_dgeikmm", (DL_FUNC) &_telefit_r_dgeikmm, 3},
    {"r_ll",                        (DL_FUNC) &r_ll,                        20},
    {"r_stpcomposition",            (DL_FUNC) &r_stpcomposition,            27},
    {"r_stpfit",                    (DL_FUNC) &r_stpfit,                    36},
    {"r_svcfit",                    (DL_FUNC) &r_svcfit,                    18},
    {"r_svcpredict",                (DL_FUNC) &r_svcpredict,                10},
    {NULL, NULL, 0}
};

RcppExport void R_init_telefit(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
