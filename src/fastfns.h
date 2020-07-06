#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Make these functions available *within* C++ in addition to R
double fastMean(const arma::vec & X);
double fastVar(const arma::vec & X);
arma::vec fastcolMeans(const arma::mat & X);
arma::vec fastrowMeans(const arma::mat & X);
arma::vec fastrowMads(const arma::mat & X);
arma::vec fastrowMads2(const arma::mat & X);
arma::vec fastrowMaxs(const arma::mat & X);
SEXP fastrowMaxs2(SEXP x);
arma::mat fastcova(const arma::mat & X);
Rcpp::NumericVector zscore_vec(const arma::vec & X);
arma::mat vecdist_arma(const arma::vec & x);
arma::mat vecdist3(const arma::vec & A);
Rcpp::List fast_dists_and_trans(const arma::vec & x, const arma::vec & Z);
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector & x);
void minus_c(double f[],double &x,double *y,int offset,int &len);
Rcpp::NumericVector avg_rank(const Rcpp::NumericVector & x);
Rcpp::List fast_dists_and_trans_by_unit(const Rcpp::NumericVector & x, const Rcpp::NumericVector & Z);
double fastmad(const Rcpp::NumericVector & x);
void inner_dist(double f[],double &x,double *y,int offset,int &len);
Rcpp::List fast_dists_and_trans_by_unit_arma(const arma::vec & x,const arma::vec & Z);
Rcpp::List fast_dists_and_trans_by_unit_arma2(const arma::vec & x,const arma::vec & Z);
double fastmad_arma(const arma::vec & x);
Rcpp::List fast_dists_by_unit_arma2_par(const arma::vec & x,const arma::vec & Z, const int threads);

