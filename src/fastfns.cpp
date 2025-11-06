// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppDist.h>
#include <cmath>
#include <iostream>
#include <vector>

#include "fastfns.h"

// #include <algorithm>
// #include <iterator>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace arma;
using namespace Rcpp;

extern "C" {
SEXP vecdist(SEXP x);
}












//[[Rcpp::export]]
arma::mat vecdist_arma(const arma::vec &x) {
  int n = x.n_elem, i = 0, j = 0;
  double d;
  arma::mat out(n, n);
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      // d = sqrt(pow( (x(i)-x(j)),2.0));
      d = std::abs(x(i) - x(j));
      out(j, i) = d;
      out(i, j) = d;
    }
  }
  return (out);
}




//[[Rcpp::export]]
arma::mat vecdist3_arma(const arma::vec &A) {
  /// This is just to compute the pairwise distances of a vector
  // Z is not used here. But in the API for distance functions.
  //  Calc euclidean distance of x -- return a square matrix n x n with all
  //  distances
  //  http://nonconditional.com/2014/04/on-the-trick-for-computing-the-squared-euclidian-distances-between-two-s
  //  ets-of-vectors/ see also
  //  https://stackoverflow.com/questions/35273292/eigen-calculate-the-distance-matrix-between-2-sets-o
  //  f-vectors and
  //  https://github.com/RoyiAvital/Projects/tree/master/CalcDistanceMatrix
  // int n = A.nrow();
  // arma::mat A = arma::mat(x.begin(), n, 1, false);
  arma::uword n = A.n_elem;

  // arma::colvec An =  sum(square(A),1);
  arma::mat C(n, n);
  arma::colvec An(n);

  An = square(A);
  C = -2 * (A * A.t());
  C.each_col() += An;
  C.each_row() += An.t();

  return sqrt(C);
}

// from mn.h in Rfast void minus_c(double f[],double &,double *,int,int &);
// and from mn.cpp
// vecdist
void minus_c(double f[], double &x, double *y, int offset, int &len) {
  double *ff = f;
  for (int i = 0; i < len; ++i, ff += offset, ++y) {
    *ff = std::abs(x - *y);
  }
}

// from vecdist.cpp in Rfast adapted here
//[[Rcpp::export]]
Rcpp::NumericMatrix vecdist2(const Rcpp::NumericVector &x) {
  int nrow, ncol;
  nrow = ncol = x.length();
  SEXP F = PROTECT(Rf_allocMatrix(REALSXP, ncol, nrow));
  double *xx = REAL(x), *end = xx + nrow, *f = REAL(F), *y = xx;
  for (; xx != end; ++xx, f += nrow)
    minus_c(f, *xx, y, 1, ncol);
  UNPROTECT(1);
  return F;
}


/*  Average (mid-) ranks like R's rank(..., ties.method = "average").
 *  Returns a double vector of length n, 1-based.
 *
 *  No special NA handling: if x contains NaN you will get NaN ranks.
 */
// [[Rcpp::export]]
arma::vec avg_rank_arma(const arma::vec &x) {

  const arma::uword n = x.n_elem;
  if (n == 0) return arma::vec();          // nothing to do

  // order vector: positions 0..n-1 sorted by ascending x
  arma::uvec ord = sort_index(x, "ascend");

  arma::vec rankv(n);

  for (arma::uword i = 0; i < n; ) {
    // find end of tie block
    arma::uword j = i + 1;
    while (j < n && x[ord[i]] == x[ord[j]]) ++j;

    double mid = (i + j - 1) / 2.0 + 1.0;   // average rank, 1-based
    for (arma::uword k = i; k < j; ++k)
      rankv[ord[k]] = mid;

    i = j;
  }
  return rankv;
}



//arma::vec avg_rank_arma(const arma::vec &x) {
//  R_xlen_t sz = x.n_elem;
//  arma::vec w = arma::linspace(0, sz - 1, sz);
//  std::sort(w.begin(), w.end(),
//            [&x](std::size_t i, std::size_t j) { return x[i] < x[j]; });
//
//  arma::vec r(sz);
//  for (R_xlen_t n, i = 0; i < sz; i += n) {
//    n = 1;
//    while (i + n < sz && x[w[i]] == x[w[i + n]])
//      ++n;
//    for (R_xlen_t k = 0; k < n; k++) {
//      r[w[i + k]] = i + (n + 1) / 2.;
//    }
//  }
//  return r;
//}

// see http://arma.sourceforge.net/docs.html#conv_to
typedef std::vector<double> stdvec;

template <typename T> Rcpp::NumericVector arma2vec(const T &x) {
  // https://stackoverflow.com/questions/14253069/convert-rcpparmadillo-vector-to-rcpp-vector
  return Rcpp::NumericVector(x.begin(), x.end());
}
// [[Rcpp::export]]
Rcpp::NumericVector avg_rank(const Rcpp::NumericVector &x) {
  R_xlen_t sz = x.size();
  Rcpp::IntegerVector w = Rcpp::seq(0, sz - 1);
  // Don't need to worry about NAs so commenting this out
  // std::sort(w.begin(), w.end(), Comparator(x));
  std::sort(w.begin(), w.end(),
            [&x](std::size_t i, std::size_t j) { return x[i] < x[j]; });

  Rcpp::NumericVector r = Rcpp::no_init_vector(sz);
  for (R_xlen_t n, i = 0; i < sz; i += n) {
    n = 1;
    while (i + n < sz && x[w[i]] == x[w[i + n]])
      ++n;
    for (R_xlen_t k = 0; k < n; k++) {
      r[w[i + k]] = i + (n + 1) / 2.;
    }
  }

  return r;
}

// [[Rcpp::export]]
double fastmad_arma(const arma::vec &x) {
  // https://gallery.rcpp.org/articles/robust-estimators/
  // scale_factor = 1.4826; default for normal distribution consistent with R
  double medx = arma::median(x);
  arma::vec med_devs = arma::abs(x - medx);
  return (1.4826 * arma::median(med_devs));
}

// from edited from Rfast dists_vec.cpp


// these next from
// https://stackoverflow.com/questions/74818733/translate-outer-from-base-r-to-rcpparmadillo




// [[Rcpp::export]]
double fastmad(const Rcpp::NumericVector &x, double center) {
  // Calculate the median of the vector if the center is not specified.

  if (center < 0) {
    center = median(x);
  }

  // Calculate the median absolute distance from the center.
  NumericVector med_devs = Rcpp::abs(x - center);
  return (1.4826 * median(med_devs));
}


