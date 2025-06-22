// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Make these functions available *within* C++ in addition to R
double fastMean(const arma::vec &X);
double fastVar(const arma::vec &X);
arma::vec fastcolMeans(const arma::mat &X);
arma::vec fastrowMeans(const arma::mat &X);
arma::vec fastrowMads(const arma::mat &X);
arma::vec fastrowMads2(const arma::mat &X);
arma::vec fastrowMaxs(const arma::mat &X);
SEXP fastrowMaxs2(SEXP x);
arma::mat fastcova(const arma::mat &X);
Rcpp::NumericVector zscore_vec(const arma::vec &X);
arma::mat vecdist_arma(const arma::vec &x);
arma::mat vecdist_squared(const arma::vec &A);
arma::mat vecdist3_arma(const arma::vec &A);
Rcpp::List fast_dists_and_trans(const arma::vec &x);
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector &x);
void minus_c(double f[], double &x, double *y, int offset, int &len);
Rcpp::NumericVector avg_rank(const Rcpp::NumericVector &x);
arma::vec avg_rank_arma(const arma::vec &x);
Rcpp::List fast_dists_and_trans_by_unit(const Rcpp::NumericVector &x);
double fastmad(const Rcpp::NumericVector &x, double center = -99.0);
void inner_dist(double f[], double &x, double *y, int offset, int &len);
Rcpp::List fast_dists_and_trans_by_unit_arma(const arma::vec &x);
Rcpp::List fast_dists_and_trans_by_unit_arma2(const arma::vec &x);
double fastmad_arma(const arma::vec &x);
Rcpp::List fast_dists_and_trans_by_unit_arma2_par(const arma::vec &x, const int threads);
Rcpp::List fast_dists_and_trans_new(const arma::vec &x);
Rcpp::List fast_dists_and_trans_new_omp(const arma::vec &x, const int threads);
Rcpp::List fast_dists_and_trans_sort(const arma::vec &x);
Rcpp::List fast_dists_and_trans_hybrid(const arma::vec &x);
Rcpp::List fast_dists_and_trans_nomax_hybrid(const arma::vec &x);
