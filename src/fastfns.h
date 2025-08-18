// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// Make these functions available *within* C++ in addition to R
arma::mat vecdist_arma(const arma::vec &x);
arma::mat vecdist3_arma(const arma::vec &A);
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector &x);
void minus_c(double f[], double &x, double *y, int offset, int &len);
Rcpp::NumericVector avg_rank(const Rcpp::NumericVector &x);
arma::vec avg_rank_arma(const arma::vec &x);
double fastmad(const Rcpp::NumericVector &x, double center = -99.0);
double fastmad_arma(const arma::vec &x);
Rcpp::List fast_dists_and_trans_hybrid(const arma::vec &x);
Rcpp::List fast_dists_and_trans_nomax_hybrid(const arma::vec &x);
