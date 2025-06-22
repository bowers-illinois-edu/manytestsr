// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include "fastfns.h"

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans(const arma::vec &x) {
  // https://stackoverflow.com/questions/24997510/how-to-speed-up-this-rcpp-function
  //  Z is not used
  int n = x.n_elem;
  arma::mat dx(n, n), dxRank0(n, n);
  arma::vec mndx(n), mndrx(n), maxdx(n), rankx(n), tanhx(n);

  rankx = avg_rank_arma(x);

  dx = vecdist3_arma(x);
  dxRank0 = vecdist3_arma(rankx);
  // Rcpp::NumericMatrix drxTmp =
  // vecdist2(Rcpp::NumericVector(rankx.begin(),rankx.end())); arma::mat
  // dxRank0(drxTmp.begin(), n, n, false);

  mndx = mean(dx, 1);
  mndrx = mean(dxRank0, 1);
  // maddx = fastrowMads2(dx);
  // maddrx = fastrowMads2(dxRank0);
  maxdx = fastrowMaxs(dx);
  // maxdrx = fastrowMaxs(dxRank0);
  // zx = zscore_vec(x);
  tanhx = tanh(x);
  // Order here matters:
  // outcome_names <-c(theresponse, "mndist", "mndistRank0", "maxdist", "rankY",
  // "tanhx")
  List res = List::create(mndx, mndrx, maxdx, rankx, tanhx);
  // List L = List::create(Named("name1") = v1 , _["name2"] = v2);
  return (res);
}

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans_by_unit_arma(const arma::vec &x) {
  // https://hydroecology.net/using-r-and-cpp-together/
  // In the end, Rcpp stuff is as efficient as using pointers etc
  arma::uword n = x.n_elem, i = 0, j = 0;

  // Setup pointers to the vector objects to make the loops faster?
  // http://adv-r.had.co.nz/C-interface.html#c-vectors
  // arma::vec mndx = Rcpp::no_init_vector(n);
  arma::vec mndx(n);
  arma::vec maxdx(n);
  // arma::vec maddx(n);
  arma::vec mndrx(n);
  // arma::vec maddrx(n);
  arma::vec maxdrx(n);
  arma::vec dist_i(n);
  // arma::vec dist_i(distptr,n,1,false,true);
  arma::vec dist_rank_i(n);
  // arma::vec zx(n);
  arma::vec rankx(n);
  arma::vec tanhx(n);

  // zx = zscore_vec(x);
  // zx.elem(find_nonfinite(zx)).zeros();

  rankx = avg_rank_arma(x);
  tanhx = tanh(x);

  // Make the distances between each unit and each other unit and their
  // summaries But don't use a whole matrix since we don't need it.
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      dist_i(j) = std::abs(x(i) - x(j));
      dist_rank_i(j) = std::abs(rankx(i) - rankx(j));
    }
    // Rcpp::Rcout << "dist_i" << std::endl;
    // Rcpp::Rcout << dist_i << std::endl;

    mndx(i) = arma::mean(dist_i);
    maxdx(i) = arma::max(dist_i);
    mndrx(i) = arma::mean(dist_rank_i);
    // maxdrx(i) = arma::max(dist_rank_i);
    // maddx(i) = fastmad_arma(dist_i);
    // maddrx(i) = fastmad_arma(dist_rank_i);
  }

  // Order here matters:
  // outcome_names <-c(theresponse, "mndist", "mndistRank0", "maxdist", "rankY",
  // "tanhx")
  Rcpp::List res = List::create(mndx, mndrx, maxdx, rankx, tanhx);
  return (res);
}

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans_by_unit_arma2_par(const arma::vec &x,
                                                  const int threads = 4) {
#ifdef _OPENMP
  if (threads > 0)
    omp_set_num_threads(threads);
    // REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif

  arma::uword n = x.n_elem, i;

  // Setup pointers to the vector objects to make the loops faster?
  // http://adv-r.had.co.nz/C-interface.html#c-vectors
  // arma::vec mndx = Rcpp::no_init_vector(n);
  arma::vec mndx(n);
  arma::vec maxdx(n);
  arma::vec mndrx(n);
  // arma::vec maxdrx(n);
  // arma::vec maddx(n);
  // arma::vec maddrx(n);
  arma::vec dist_i(n);
  arma::vec dist_rank_i(n);

  // arma::vec zx(n);
  arma::vec rankx(n);
  arma::vec tanhx(n);

  // zx = zscore_vec(x);
  // zx.elem(find_nonfinite(zx)).zeros();

  rankx = avg_rank_arma(x);
  tanhx = tanh(x);

#pragma omp parallel for default(shared)                                       \
    schedule(static) private(i, dist_i, dist_rank_i)
  for (i = 0; i < n; i++) {
    dist_i = abs(x(i) - x);
    // dist_i = sqrt(square( x(i) - x));
    // dist_i = square( x(i) - x);
    dist_rank_i = abs(rankx(i) - rankx);
    // dist_rank_i = sqrt(square( rankx(i) - rankx));
    // dist_rank_i = square( rankx(i) - rankx);
    mndx(i) = arma::mean(dist_i);
    maxdx(i) = arma::max(dist_i);
    mndrx(i) = arma::mean(dist_rank_i);
    // maxdrx(i) = arma::max(dist_rank_i);
    // maddx(i) = fastmad_arma(dist_i);
    // maddrx(i) = fastmad_arma(dist_rank_i);
  }

  // Order here matters:
  // outcome_names <-c(theresponse, "mndist", "mndistRank0", "maxdist", "rankY",
  // "tanhx")
  Rcpp::List res = List::create(mndx, mndrx, maxdx, rankx, tanhx);
  return (res);
}

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans_by_unit_arma2(const arma::vec &x) {
  // this version avoids the inner loop using vectors.
  arma::uword n = x.n_elem, i;

  // arma::vec zx(n);
  arma::vec rankx(n);
  arma::vec tanhx(n);

  // zx = zscore_vec(x);
  // zx.elem(find_nonfinite(zx)).zeros();

  rankx = avg_rank_arma(x);
  tanhx = tanh(x);

  // Setup pointers to the vector objects to make the loops faster?
  // http://adv-r.had.co.nz/C-interface.html#c-vectors
  // arma::vec mndx = Rcpp::no_init_vector(n);
  arma::vec mndx(n);
  arma::vec maxdx(n);
  arma::vec mndrx(n);
  arma::vec maxdrx(n);
  arma::vec maddx(n);
  arma::vec maddrx(n);

  for (i = 0; i < n; i++) {
    arma::vec dist_i(n);
    arma::vec dist_rank_i(n);
    // dist_i = sqrt(pow( x(i) - x,2 ));
    // dist_rank_i = sqrt(pow( rankx(i) - rankx, 2 ));
    // dist_i = pow( x(i) - x,2 );
    dist_i = arma::abs(x(i) - x);
    // dist_rank_i = arma::pow( rankx(i) - rankx, 2 );
    dist_rank_i = arma::abs(rankx(i) - rankx);
    mndx(i) = arma::mean(dist_i);
    maxdx(i) = arma::max(dist_i);
    // maddx(i) = fastmad_arma(dist_i);
    mndrx(i) = arma::mean(dist_rank_i);
    // maxdrx(i) = arma::max(dist_rank_i);
    // maddrx(i) = fastmad_arma(dist_rank_i);
  }
  // Rcpp::Rcout << "dist_i" << std::endl;
  // Rcpp::Rcout << dist_i << std::endl;

  // Order here matters:
  // outcome_names <-c(theresponse, "mndist", "mndistRank0", "maxdist", "rankY",
  // "tanhx")
  Rcpp::List res = List::create(mndx, mndrx, maxdx, rankx, tanhx);
  return (res);
}


/*  Compute, for every unit i:
 *    – mean |x_i − x_j|               (col 1)
 *    – mean |rank_i − rank_j|         (col 2)
 *    – max  |x_i − x_j|               (col 3)
 *    – raw rank (average-rank)        (col 4)
 *    – tanh(x_i)                      (col 5)
 *
 *  All distances are taken across *all other* units j ≠ i.
 *
 *  @param x Numeric vector of outcomes.
 *  @param Z Unused at the moment (reserved for future stratification).
 *  @return R list with named components.
 */
// [[Rcpp::export]]
Rcpp::List fast_dists_and_trans_new(const arma::vec &x) {

  const arma::uword n = x.n_elem;
  if (n < 2) stop("x must have length >= 2");

  arma::vec rankx = avg_rank_arma(x);      // mid-ranks, ties OK
  arma::vec tanhx = arma::tanh(x);

  arma::vec mean_raw(n,  arma::fill::zeros);
  arma::vec mean_rank(n, arma::fill::zeros);
  arma::vec max_raw(n,   arma::fill::zeros);
//  arma::vec max_rank(n,  arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    double xi = x[i], ri = rankx[i];
    double sum_raw  = 0.0, sum_rank = 0.0;
    double mx_raw   = 0.0; //, mx_rank  = 0.0;

    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;
      double d  = std::abs(xi - x[j]);
      double dr = std::abs(ri - rankx[j]);

      sum_raw  += d;
      sum_rank += dr;
      if (d  > mx_raw)  mx_raw  = d;
 //     if (dr > mx_rank) mx_rank = dr;
    }
    mean_raw[i]  = sum_raw  / (n - 1);
    mean_rank[i] = sum_rank / (n - 1);
    max_raw[i]   = mx_raw;
//    max_rank[i]  = mx_rank;
  }

  return List::create(
    _["mean_dist"]        = mean_raw,
    _["mean_rank_dist"]   = mean_rank,
    _["max_dist"]         = max_raw,
 //   _["max_rank_dist"]    = max_rank,
    _["rankY"]             = rankx,
    _["tanhY"]             = tanhx
  );
}



/*  Compute, for every unit i, using parallelism:
 *    – mean |x_i − x_j|               (col 1)
 *    – mean |rank_i − rank_j|         (col 2)
 *    – max  |x_i − x_j|               (col 3)
 *    – raw rank (average-rank)        (col 4)
 *    – tanh(x_i)                      (col 5)
 *
 *  All distances are taken across *all other* units j ≠ i.
 *
 *  @param x Numeric vector of outcomes.
 *  @param Z Unused at the moment (reserved for future stratification).
 *  @return R list with named components.
 */
// [[Rcpp::export]]
Rcpp::List fast_dists_and_trans_new_omp(const arma::vec &x, const int threads=0) {

#ifdef _OPENMP
  // turn off OpenMP’s own dynamic adjustment (optional but predictable)
  omp_set_dynamic(0);

  if (threads > 0)
    threads = std::min<int>(threads, omp_get_max_threads());
  omp_set_num_threads(threads);            // user-supplied
  /* else              */
  /*   keep whatever OMP_NUM_THREADS or omp_set_num_threads()  */
  /*   has already established for this R session              */
#endif

  const arma::uword n = x.n_elem;
  if (n < 2) stop("x must have length >= 2");

  arma::vec rankx = avg_rank_arma(x);      // mid-ranks, ties OK
  arma::vec tanhx = arma::tanh(x);

  arma::vec mean_raw(n,  arma::fill::zeros);
  arma::vec mean_rank(n, arma::fill::zeros);
  arma::vec max_raw(n,   arma::fill::zeros);
  //  arma::vec max_rank(n,  arma::fill::zeros);

#ifdef _OPENMP
  #pragma omp parallel for default(none) shared(n,x,rankx,mean_raw,mean_rank,max_raw,max_rank)
#endif

  for (arma::uword i = 0; i < n; ++i) {
    double xi = x[i], ri = rankx[i];
    double sum_raw  = 0.0, sum_rank = 0.0;
    double mx_raw   = 0.0; //, mx_rank  = 0.0;

    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;
      double d  = std::abs(xi - x[j]);
      double dr = std::abs(ri - rankx[j]);

      sum_raw  += d;
      sum_rank += dr;
      if (d  > mx_raw)  mx_raw  = d;
      //    if (dr > mx_rank) mx_rank = dr;
    }
    mean_raw[i]  = sum_raw  / (n - 1);
    mean_rank[i] = sum_rank / (n - 1);
    max_raw[i]   = mx_raw;
    // max_rank[i]  = mx_rank;
  }

  return List::create(
    _["mean_dist"]        = mean_raw,
    _["mean_rank_dist"]   = mean_rank,
    _["max_dist"]         = max_raw,
  //  _["max_rank_dist"]    = max_rank,
    _["rankY"]             = rankx,
    _["tanhY"]             = tanhx
  );
}


//' Fast per-unit distance summaries (scalar outcome)
//'
//' Computes, for **each observation** \eqn{i = 1,\dots,n} in a numeric vector
//' \eqn{x}, scalar summaries plus simple transforms, **without ever
//' forming the full \eqn{n\times n} distance matrix**.
//'
//' \itemize{
//'   \item **`mean_dist`** – mean absolute distance
//'         \deqn{\frac{1}{n-1}\sum_{j\neq i}|x_i-x_j|}{}
//'   \item **`mean_rank_dist`** – same mean on the mid-rank scale;
//'         closed-form, no second loop.
//'   \item **`max_dist`** – maximum absolute distance
//'         \eqn{\max\{\,x_i-\min(x),\;\max(x)-x_i\,\}}{}
//'   \item **`rankY`** – average (mid-) rank of \code{x} (\code{ties="average"}).
//'   \item **`tanhY`** – element-wise \eqn{\tanh(x_i)} shrink transform.
//' }
//'
//' **Complexity**
//' \itemize{
//'   \item \eqn{O(n \log n)} time (sorting + prefix sums)
//'   \item \eqn{O(n)}  space (only vectors of length \code{n})
//' }
//' For \code{n = 10\,000} this typically runs in ≈2 ms on an Apple M-series
//' core with < 0.5 MB peak RAM, much faster and far lighter than allocating
//' the full distance matrix.
//'
//' @param x Numeric vector – the scalar outcome for all units.
//'
//' @return A named \code{list} with components
//'   \describe{
//'     \item{\code{mean_dist}}{numeric vector, length \code{n}.}
//'     \item{\code{mean_rank_dist}}{numeric vector, length \code{n}.}
//'     \item{\code{max_dist}}{numeric vector, length \code{n}.}
//'     \item{\code{rankY}}{average ranks (mid-ranks).}
//'     \item{\code{tanhY}}{\eqn{\tanh(x_i)} values.}
//'   }
//'
//' @examples
//' set.seed(1)
//' x <- rnorm(8)
//' fast_dists_and_trans_hybrid(x)
//'
//' ## compare to explicit distance matrix (slow / big):
//' dx <- abs(outer(x, x, "-"))
//' mean_dist_ref <- colSums(dx) / (length(x) - 1)
//' stopifnot(all.equal(fast_dists_and_trans_hybrid(x)$mean_dist,
//'                     mean_dist_ref))
//'
//' @name fast_dists_and_trans_hybrid
//' @export
// [[Rcpp::export]]
Rcpp::List fast_dists_and_trans_hybrid(const arma::vec &x) {

  const arma::uword n = x.n_elem;
  if (n < 2) stop("x must have length >= 2");

  // ----- ranks (original order) -------------------------------------------
  arma::vec rankx = avg_rank_arma(x);               // mid-ranks, ties OK

  // ----- sort x once + prefix sums  ---------------------------------------
  arma::uvec ord = sort_index(x);                   // indices that sort x
  arma::vec  xs  = x(ord);                          // sorted copy
  arma::vec  ps  = cumsum(xs);                      // prefix sums

  // ----- mean raw distances (O(n log n)) ----------------------------------
  arma::vec mean_raw_sorted(n);
  for (arma::uword k = 0; k < n; ++k) {
    double left  = xs[k] * k - (k ? ps[k-1] : 0.0);
    double right = (ps[n-1] - ps[k]) - xs[k] * (n - k - 1);
    mean_raw_sorted[k] = (left + right) / (n - 1);
  }
  arma::vec mean_raw(n);
  mean_raw.elem(ord) = mean_raw_sorted;       // back to data order

  // ----- mean rank distances (handles ties) ---------------------------------
  arma::uvec ord_r = sort_index(rankx);          // indices that sort the rank vector
  arma::vec  rs    = rankx(ord_r);               // sorted ranks (with duplicates)
  arma::vec  ps_r  = cumsum(rs);                 // prefix sums of ranks

  arma::vec mean_rank_sorted(n);
  for (arma::uword k = 0; k < n; ++k) {
    double left  = rs[k] * k - (k ? ps_r[k - 1] : 0.0);
    double right = (ps_r[n - 1] - ps_r[k]) - rs[k] * (n - k - 1);
    mean_rank_sorted[k] = (left + right) / (n - 1);
  }
  arma::vec mean_rank(n);
  mean_rank.elem(ord_r) = mean_rank_sorted;   // back to data order

// this next requires no ties and all unique ranks
//  ----- mean rank distances (closed form) --------------------------------
//  arma::vec mean_rank(n);
//  for (arma::uword i = 0; i < n; ++i) {
//    double r = rankx[i];
//    mean_rank[i] = ((r - 1.0) * r + (n - r) * (n - r + 1.0)) / 2.0;
//    mean_rank[i] /= (n - 1);
//  }

  // ----- max raw & rank distances (O(n)) ----------------------------------
  double xmin = xs[0], xmax = xs[n-1];
  arma::vec max_raw(n);
  for (arma::uword i = 0; i < n; ++i) {
    double xi = x[i];
    max_raw[i]  = std::max(xi - xmin, xmax - xi);
  }

  arma::vec tanhx = arma::tanh(x);   // column vector, length n
  // ----- wrap and return ---------------------------------------------------
  return List::create(
    _["mean_dist"]        = mean_raw,
    _["mean_rank_dist"]   = mean_rank,
    _["max_dist"]         = max_raw,
    _["rankY"]            = rankx,
    _["tanhY"]            = tanhx
  );
}



// [[Rcpp::export]]
Rcpp::List fast_dists_and_trans_nomax_hybrid(const arma::vec &x) {

  const arma::uword n = x.n_elem;
  if (n < 2) stop("x must have length >= 2");

  // ----- ranks (original order) -------------------------------------------
  arma::vec rankx = avg_rank_arma(x);               // mid-ranks, ties OK

  // ----- sort x once + prefix sums  ---------------------------------------
  arma::uvec ord = sort_index(x);                   // indices that sort x
  arma::vec  xs  = x(ord);                          // sorted copy
  arma::vec  ps  = cumsum(xs);                      // prefix sums

  // ----- mean raw distances (O(n log n)) ----------------------------------
  arma::vec mean_raw_sorted(n);
  for (arma::uword k = 0; k < n; ++k) {
    double left  = xs[k] * k - (k ? ps[k-1] : 0.0);
    double right = (ps[n-1] - ps[k]) - xs[k] * (n - k - 1);
    mean_raw_sorted[k] = (left + right) / (n - 1);
  }
  arma::vec mean_raw(n);
  mean_raw.elem(ord) = mean_raw_sorted;       // back to data order

  // ----- mean rank distances (closed form) --------------------------------
  arma::vec mean_rank(n);
  for (arma::uword i = 0; i < n; ++i) {
    double r = rankx[i];
    mean_rank[i] = ((r - 1.0) * r + (n - r) * (n - r + 1.0)) / 2.0;
    mean_rank[i] /= (n - 1);
  }

  arma::vec tanhx = arma::tanh(x);   // column vector, length n
  // ----- wrap and return ---------------------------------------------------
  return List::create(
    _["mean_dist"]        = mean_raw,
    _["mean_rank_dist"]   = mean_rank,
    _["rankY"]            = rankx,
    _["tanhY"]            = tanhx
  );
}

