// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

#include "fastfns.h"

//[[Rcpp::export]]
Rcpp::List fast_dists_and_trans(const arma::vec &x, const arma::vec &Z) {
  // https://stackoverflow.com/questions/24997510/how-to-speed-up-this-rcpp-function
  //  Z is not used
  int n = x.n_elem;
  arma::mat dx(n, n), dxRank0(n, n);
  arma::vec mndx(n), mndrx(n), maxdx(n), rankx(n), tanhx(n);

  dx = vecdist3_arma(x);
  //  Rcpp::NumericMatrix dxTmp =
  //  vecdist2(Rcpp::NumericVector(x.begin(),x.end()));
  // dx = as<arma::mat>(dxTmp);
  // arma::mat dx(dxTmp.begin(), n, n, false);

  rankx = avg_rank_arma(x);
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
Rcpp::List fast_dists_and_trans_by_unit_arma(const arma::vec &x,
                                             const arma::vec &Z) {
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
                                                  const arma::vec &Z,
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
Rcpp::List fast_dists_and_trans_by_unit_arma2(const arma::vec &x,
                                              const arma::vec &Z) {
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
