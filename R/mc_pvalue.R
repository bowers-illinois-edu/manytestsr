## Monte Carlo p-value convention for the small-sample permutation branch.
##
## Below simthresh the p-value functions in pval_fns.R draw their reference
## distribution with coin::approximate(nresample = sims). coin reports
## count/sims, where count is the number of resampled statistics at least as
## extreme as the observed one. When no resample beats the observed value that
## is exactly zero -- reproduced directly with perfect separation at
## sims = 9 and sims = 49.
##
## A p-value of zero sits at or below every threshold, so the node always
## rejects. Under the null the observed statistic is the most extreme of the
## sims + 1 exchangeable values with probability about 1/(sims + 1), so any
## critical value below 1/(sims + 1) has realized size about 1/(sims + 1)
## rather than its nominal value. The damage is worst for the smallest
## thresholds, which is to say for alpha/m-style corrections -- so leaving it
## in place would bias a head-to-head comparison in favour of procedures that
## test near nominal alpha, this package's included.
##
## The observed assignment is itself a draw from the null and belongs in both
## the numerator and the denominator. Reporting (count + 1)/(sims + 1) is the
## standard convention (Davison and Hinkley 1997, Sec. 4.2; Phipson and Smyth
## 2010, "Permutation P-values should never be zero"), is never
## anti-conservative, and cannot return zero.

#' Apply the Monte Carlo p-value convention
#'
#' Converts the count/sims p-value that \code{coin} reports on the
#' simulation branch into \eqn{(count + 1) / (sims + 1)}. Asymptotic
#' p-values are returned unchanged.
#'
#' @param p Numeric p-value as returned by \code{coin::pvalue()}.
#' @param sims Number of resamples passed to \code{coin::approximate()}.
#' @param thedist The distribution argument handed to the \code{coin} test:
#'   either the string \code{"asymptotic"} or an approximate null
#'   distribution. Used to decide whether the correction applies.
#' @return The corrected p-value, in \eqn{[1/(sims + 1), 1]} on the
#'   simulation branch and unchanged on the asymptotic branch.
#' @keywords internal
adjust_mc_pvalue <- function(p, sims, thedist) {
  ## The asymptotic branch has no resampling and needs no correction.
  if (identical(thedist, "asymptotic")) {
    return(p)
  }
  if (is.na(p)) {
    return(p)
  }
  ## coin reports count/sims exactly, so the count recovers cleanly. round()
  ## absorbs the floating-point error in that division.
  count <- round(p * sims)
  return((count + 1) / (sims + 1))
}
