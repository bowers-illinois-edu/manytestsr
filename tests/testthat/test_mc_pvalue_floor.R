# Tests for the Monte Carlo p-value convention on the small-sample branch.
# Written 2026-08-14.
#
# THE DEFECT. Below simthresh the reference distribution is Monte Carlo
# permutation via coin::approximate(nresample = sims). coin returns
# count/sims, where count is the number of resampled statistics at least as
# extreme as the observed one. When no resample beats the observed value that
# is exactly 0 -- verified directly: with perfect separation coin returns
# p = 0.000000 at sims = 9 and at sims = 49.
#
# WHY IT MATTERS HERE. A p-value of 0 is at or below every threshold, so such
# a node always rejects. Under the null the observed statistic is the most
# extreme of the sims + 1 exchangeable values with probability about
# 1/(sims + 1), so any critical value below 1/(sims + 1) has realized size
# about 1/(sims + 1) rather than its nominal value. That is anti-conservative
# exactly at the small nodes where this branch runs, and exactly for the
# small thresholds that alpha/m-style corrections need -- so it biases
# COMPETITOR procedures more than the top-down one, which tests near nominal
# alpha. It must be fixed before any head-to-head benchmark.
#
# THE FIX. Report (count + 1)/(sims + 1) instead of count/sims: the observed
# assignment is itself a draw from the null and belongs in both the numerator
# and the denominator. This is the standard convention (Davison and Hinkley
# 1997; Phipson and Smyth 2010, "Permutation P-values should never be zero").
# It is never anti-conservative and it can never return 0.
#
# These tests FAIL before the fix and pass after it.

test_that("the MC branch never returns a p-value of exactly zero", {
  # Perfect separation: the observed statistic is the maximum of the
  # permutation distribution, so with few resamples none can beat it. This is
  # the configuration that produced p = 0 from coin directly.
  sepdat <- data.frame(
    blockF = factor(rep("b1", 12)),
    trtF = factor(rep(c(0L, 1L), each = 6)),
    Y1 = c(1:6, 101:106)
  )
  for (B in c(9, 19, 49)) {
    p_ow <- pOneway(
      fmla = Y1 ~ trtF | blockF, dat = sepdat,
      simthresh = 20, sims = B, parallel = "no"
    )
    p_wc <- pWilcox(
      fmla = Y1 ~ trtF | blockF, dat = sepdat,
      simthresh = 20, sims = B, parallel = "no"
    )
    expect_gt(p_ow, 0)
    expect_gt(p_wc, 0)
    # The smallest attainable value is exactly 1/(B + 1).
    expect_gte(p_ow, 1 / (B + 1) - 1e-12)
    expect_gte(p_wc, 1 / (B + 1) - 1e-12)
  }
})

test_that("MC p-values land on the (count + 1)/(sims + 1) lattice", {
  set.seed(20260814)
  B <- 99
  # Any MC p-value must be k/(B + 1) for an integer k in 1..(B + 1).
  for (i in 1:15) {
    d <- data.frame(
      blockF = factor(rep("b1", 12)),
      trtF = factor(rep(c(0L, 1L), each = 6)),
      Y1 = rnorm(12)
    )
    p <- pOneway(
      fmla = Y1 ~ trtF | blockF, dat = d,
      simthresh = 20, sims = B, parallel = "no"
    )
    k <- p * (B + 1)
    expect_lt(abs(k - round(k)), 1e-8)
    expect_gte(round(k), 1)
    expect_lte(round(k), B + 1)
  }
})

test_that("the MC branch is valid at its own smallest threshold", {
  # Under the sharp null, P(p <= t) <= t must hold at t = 1/(sims + 1), the
  # threshold where the old convention was worst. With count/sims the realized
  # size was about 2/(sims + 1) there; with the fix it is at most 1/(sims + 1).
  set.seed(20260814)
  B <- 19
  nsim <- 2000
  t_small <- 1 / (B + 1)
  ps <- replicate(nsim, {
    d <- data.frame(
      blockF = factor(rep("b1", 10)),
      trtF = factor(rep(c(0L, 1L), each = 5)),
      Y1 = rnorm(10)
    )
    pOneway(
      fmla = Y1 ~ trtF | blockF, dat = d,
      simthresh = 20, sims = B, parallel = "no"
    )
  })
  realized <- mean(ps <= t_small + 1e-12)
  # 3 Monte Carlo standard errors of slack, one-sided: we only object to
  # anti-conservatism.
  mcse <- sqrt(t_small * (1 - t_small) / nsim)
  expect_lte(realized, t_small + 3 * mcse)
})

test_that("the asymptotic branch is untouched by the fix", {
  # Above simthresh the distribution is Strasser-Weber asymptotic and no
  # Monte Carlo correction applies. Values must be unchanged.
  set.seed(20260814)
  bigdat <- data.frame(
    blockF = factor(rep("b1", 60)),
    trtF = factor(rep(c(0L, 1L), each = 30)),
    Y1 = rnorm(60)
  )
  p_pkg <- pOneway(
    fmla = Y1 ~ trtF | blockF, dat = bigdat,
    simthresh = 20, sims = 1000, parallel = "no"
  )
  p_coin <- coin::pvalue(coin::oneway_test(
    Y1 ~ trtF | blockF,
    data = bigdat, distribution = "asymptotic"
  ))[[1]]
  expect_equal(p_pkg, as.numeric(p_coin))
})

test_that("pIndepDist also respects the convention on the MC branch", {
  # KNOWN ENVIRONMENT CRASH, not a defect in this fix. pIndepDist() calls
  # distfn() inside a data.table `:=`, and that C++ distance layer aborts R
  # on Jake's machine. Verified this session that the abort is identical on
  # the pre-fix code (git stash), so it predates the p-value convention
  # change; see HANDOFF.md on the nine test files that crash R here. The
  # correction is applied at the pIndepDist site in the same one-line form
  # as the other five, and is covered there by pOneway and pWilcox.
  skip("pIndepDist crashes R on this machine (pre-existing C++ dist layer)")
})
