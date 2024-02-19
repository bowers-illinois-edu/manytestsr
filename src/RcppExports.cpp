// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_dists_and_trans
Rcpp::List fast_dists_and_trans(const arma::vec& x, const arma::vec& Z);
RcppExport SEXP _manytestsr_fast_dists_and_trans(SEXP xSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dists_and_trans(x, Z));
    return rcpp_result_gen;
END_RCPP
}
// fast_dists_and_trans_by_unit_arma
Rcpp::List fast_dists_and_trans_by_unit_arma(const arma::vec& x, const arma::vec& Z);
RcppExport SEXP _manytestsr_fast_dists_and_trans_by_unit_arma(SEXP xSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dists_and_trans_by_unit_arma(x, Z));
    return rcpp_result_gen;
END_RCPP
}
// fast_dists_and_trans_by_unit_arma2_par
Rcpp::List fast_dists_and_trans_by_unit_arma2_par(const arma::vec& x, const arma::vec& Z, const int threads);
RcppExport SEXP _manytestsr_fast_dists_and_trans_by_unit_arma2_par(SEXP xSEXP, SEXP ZSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dists_and_trans_by_unit_arma2_par(x, Z, threads));
    return rcpp_result_gen;
END_RCPP
}
// fast_dists_and_trans_by_unit_arma2
Rcpp::List fast_dists_and_trans_by_unit_arma2(const arma::vec& x, const arma::vec& Z);
RcppExport SEXP _manytestsr_fast_dists_and_trans_by_unit_arma2(SEXP xSEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_dists_and_trans_by_unit_arma2(x, Z));
    return rcpp_result_gen;
END_RCPP
}
// fastMean
double fastMean(const arma::vec& X);
RcppExport SEXP _manytestsr_fastMean(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMean(X));
    return rcpp_result_gen;
END_RCPP
}
// fastMedian
double fastMedian(const arma::vec& X);
RcppExport SEXP _manytestsr_fastMedian(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastMedian(X));
    return rcpp_result_gen;
END_RCPP
}
// fastVar
double fastVar(const arma::vec& X);
RcppExport SEXP _manytestsr_fastVar(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastVar(X));
    return rcpp_result_gen;
END_RCPP
}
// fastcolMeans
arma::vec fastcolMeans(const arma::mat& X);
RcppExport SEXP _manytestsr_fastcolMeans(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastcolMeans(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMeans
arma::vec fastrowMeans(const arma::mat& X);
RcppExport SEXP _manytestsr_fastrowMeans(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMeans(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMads
arma::vec fastrowMads(const arma::mat& X);
RcppExport SEXP _manytestsr_fastrowMads(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMads(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMads4
Rcpp::NumericVector fastrowMads4(Rcpp::NumericMatrix& X);
RcppExport SEXP _manytestsr_fastrowMads4(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMads4(X));
    return rcpp_result_gen;
END_RCPP
}
// fastcolMads4
Rcpp::NumericVector fastcolMads4(Rcpp::NumericMatrix& X);
RcppExport SEXP _manytestsr_fastcolMads4(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastcolMads4(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMads2
arma::vec fastrowMads2(const arma::mat& X);
RcppExport SEXP _manytestsr_fastrowMads2(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMads2(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMads3
Rcpp::NumericVector fastrowMads3(const Rcpp::NumericMatrix& X);
RcppExport SEXP _manytestsr_fastrowMads3(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMads3(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMaxs
arma::vec fastrowMaxs(const arma::mat& X);
RcppExport SEXP _manytestsr_fastrowMaxs(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMaxs(X));
    return rcpp_result_gen;
END_RCPP
}
// fastrowMaxs2
SEXP fastrowMaxs2(SEXP x);
RcppExport SEXP _manytestsr_fastrowMaxs2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fastrowMaxs2(x));
    return rcpp_result_gen;
END_RCPP
}
// fastcova
arma::mat fastcova(const arma::mat& X);
RcppExport SEXP _manytestsr_fastcova(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(fastcova(X));
    return rcpp_result_gen;
END_RCPP
}
// replace_na_nan
NumericVector replace_na_nan(NumericVector x, double replacement);
RcppExport SEXP _manytestsr_replace_na_nan(SEXP xSEXP, SEXP replacementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type replacement(replacementSEXP);
    rcpp_result_gen = Rcpp::wrap(replace_na_nan(x, replacement));
    return rcpp_result_gen;
END_RCPP
}
// zscore_vec2
Rcpp::NumericVector zscore_vec2(const Rcpp::NumericVector& x);
RcppExport SEXP _manytestsr_zscore_vec2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(zscore_vec2(x));
    return rcpp_result_gen;
END_RCPP
}
// zscore_vec
Rcpp::NumericVector zscore_vec(const arma::vec& X);
RcppExport SEXP _manytestsr_zscore_vec(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(zscore_vec(X));
    return rcpp_result_gen;
END_RCPP
}
// vecdist_arma
arma::mat vecdist_arma(const arma::vec& x);
RcppExport SEXP _manytestsr_vecdist_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecdist_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// vecdist4_arma
arma::mat vecdist4_arma(const arma::vec& x);
RcppExport SEXP _manytestsr_vecdist4_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecdist4_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// vecdist3_arma
arma::mat vecdist3_arma(const arma::vec& A);
RcppExport SEXP _manytestsr_vecdist3_arma(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(vecdist3_arma(A));
    return rcpp_result_gen;
END_RCPP
}
// vecdist2
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector& x);
RcppExport SEXP _manytestsr_vecdist2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecdist2(x));
    return rcpp_result_gen;
END_RCPP
}
// avg_rank_arma
arma::vec avg_rank_arma(const arma::vec& x);
RcppExport SEXP _manytestsr_avg_rank_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(avg_rank_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// avg_rank
Rcpp::NumericVector avg_rank(const Rcpp::NumericVector& x);
RcppExport SEXP _manytestsr_avg_rank(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(avg_rank(x));
    return rcpp_result_gen;
END_RCPP
}
// fastmad_arma
double fastmad_arma(const arma::vec& x);
RcppExport SEXP _manytestsr_fastmad_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fastmad_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// manhattan_dist
NumericMatrix manhattan_dist(arma::vec x);
RcppExport SEXP _manytestsr_manhattan_dist(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(manhattan_dist(x));
    return rcpp_result_gen;
END_RCPP
}
// euc_dist_arma1
arma::mat euc_dist_arma1(const arma::vec& x);
RcppExport SEXP _manytestsr_euc_dist_arma1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(euc_dist_arma1(x));
    return rcpp_result_gen;
END_RCPP
}
// vecdist3
Rcpp::NumericMatrix vecdist3(const NumericVector& x);
RcppExport SEXP _manytestsr_vecdist3(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecdist3(x));
    return rcpp_result_gen;
END_RCPP
}
// trimmed_mean
double trimmed_mean(const NumericVector& x, double trim_percent);
RcppExport SEXP _manytestsr_trimmed_mean(SEXP xSEXP, SEXP trim_percentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type trim_percent(trim_percentSEXP);
    rcpp_result_gen = Rcpp::wrap(trimmed_mean(x, trim_percent));
    return rcpp_result_gen;
END_RCPP
}
// fastmad
double fastmad(const Rcpp::NumericVector& x, double center);
RcppExport SEXP _manytestsr_fastmad(SEXP xSEXP, SEXP centerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type center(centerSEXP);
    rcpp_result_gen = Rcpp::wrap(fastmad(x, center));
    return rcpp_result_gen;
END_RCPP
}
// fast_huberM
double fast_huberM(NumericVector x, double k, double tol, double trim);
RcppExport SEXP _manytestsr_fast_huberM(SEXP xSEXP, SEXP kSEXP, SEXP tolSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type trim(trimSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_huberM(x, k, tol, trim));
    return rcpp_result_gen;
END_RCPP
}
// col_huberM
NumericVector col_huberM(NumericMatrix dx, double k, double tol, double trim);
RcppExport SEXP _manytestsr_col_huberM(SEXP dxSEXP, SEXP kSEXP, SEXP tolSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dx(dxSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type trim(trimSEXP);
    rcpp_result_gen = Rcpp::wrap(col_huberM(dx, k, tol, trim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_manytestsr_fast_dists_and_trans", (DL_FUNC) &_manytestsr_fast_dists_and_trans, 2},
    {"_manytestsr_fast_dists_and_trans_by_unit_arma", (DL_FUNC) &_manytestsr_fast_dists_and_trans_by_unit_arma, 2},
    {"_manytestsr_fast_dists_and_trans_by_unit_arma2_par", (DL_FUNC) &_manytestsr_fast_dists_and_trans_by_unit_arma2_par, 3},
    {"_manytestsr_fast_dists_and_trans_by_unit_arma2", (DL_FUNC) &_manytestsr_fast_dists_and_trans_by_unit_arma2, 2},
    {"_manytestsr_fastMean", (DL_FUNC) &_manytestsr_fastMean, 1},
    {"_manytestsr_fastMedian", (DL_FUNC) &_manytestsr_fastMedian, 1},
    {"_manytestsr_fastVar", (DL_FUNC) &_manytestsr_fastVar, 1},
    {"_manytestsr_fastcolMeans", (DL_FUNC) &_manytestsr_fastcolMeans, 1},
    {"_manytestsr_fastrowMeans", (DL_FUNC) &_manytestsr_fastrowMeans, 1},
    {"_manytestsr_fastrowMads", (DL_FUNC) &_manytestsr_fastrowMads, 1},
    {"_manytestsr_fastrowMads4", (DL_FUNC) &_manytestsr_fastrowMads4, 1},
    {"_manytestsr_fastcolMads4", (DL_FUNC) &_manytestsr_fastcolMads4, 1},
    {"_manytestsr_fastrowMads2", (DL_FUNC) &_manytestsr_fastrowMads2, 1},
    {"_manytestsr_fastrowMads3", (DL_FUNC) &_manytestsr_fastrowMads3, 1},
    {"_manytestsr_fastrowMaxs", (DL_FUNC) &_manytestsr_fastrowMaxs, 1},
    {"_manytestsr_fastrowMaxs2", (DL_FUNC) &_manytestsr_fastrowMaxs2, 1},
    {"_manytestsr_fastcova", (DL_FUNC) &_manytestsr_fastcova, 1},
    {"_manytestsr_replace_na_nan", (DL_FUNC) &_manytestsr_replace_na_nan, 2},
    {"_manytestsr_zscore_vec2", (DL_FUNC) &_manytestsr_zscore_vec2, 1},
    {"_manytestsr_zscore_vec", (DL_FUNC) &_manytestsr_zscore_vec, 1},
    {"_manytestsr_vecdist_arma", (DL_FUNC) &_manytestsr_vecdist_arma, 1},
    {"_manytestsr_vecdist4_arma", (DL_FUNC) &_manytestsr_vecdist4_arma, 1},
    {"_manytestsr_vecdist3_arma", (DL_FUNC) &_manytestsr_vecdist3_arma, 1},
    {"_manytestsr_vecdist2", (DL_FUNC) &_manytestsr_vecdist2, 1},
    {"_manytestsr_avg_rank_arma", (DL_FUNC) &_manytestsr_avg_rank_arma, 1},
    {"_manytestsr_avg_rank", (DL_FUNC) &_manytestsr_avg_rank, 1},
    {"_manytestsr_fastmad_arma", (DL_FUNC) &_manytestsr_fastmad_arma, 1},
    {"_manytestsr_manhattan_dist", (DL_FUNC) &_manytestsr_manhattan_dist, 1},
    {"_manytestsr_euc_dist_arma1", (DL_FUNC) &_manytestsr_euc_dist_arma1, 1},
    {"_manytestsr_vecdist3", (DL_FUNC) &_manytestsr_vecdist3, 1},
    {"_manytestsr_trimmed_mean", (DL_FUNC) &_manytestsr_trimmed_mean, 2},
    {"_manytestsr_fastmad", (DL_FUNC) &_manytestsr_fastmad, 2},
    {"_manytestsr_fast_huberM", (DL_FUNC) &_manytestsr_fast_huberM, 4},
    {"_manytestsr_col_huberM", (DL_FUNC) &_manytestsr_col_huberM, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_manytestsr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
