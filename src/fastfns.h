#include <RcppArmadillo.h>

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
// arma::vec zscore_vec(const arma::vec & X);
arma::mat vecdist_arma(const arma::vec & x);
arma::mat vecdist3(const arma::vec & A);
Rcpp::List fast_dists_and_trans(const arma::vec & x, const arma::vec & Z);
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector & x);
void minus_c(double f[],double &x,double *y,int offset,int &len);

