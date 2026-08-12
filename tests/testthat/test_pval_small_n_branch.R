# Tests for the small-sample permutation branch (referee-correction Fix 4,
# see FIX_PLAN.md and FIXME_pval_permutation_branch.md). Written 2026-08-11.
#
# Context: the referee's minor concern 8 established that no coin::exact
# path exists in the package -- above simthresh the reference distribution
# is Strasser-Weber asymptotic, below it Monte Carlo permutation via
# coin::approximate(). The audit additionally found a latent defect in that
# call: a stray positional `object` argument that binds to coin's `cl`
# formal (every other formal is matched by name). Under parallel = "no"
# the promise is never forced, so the branch happens to work; under
# parallel = "snow" it would fail. The fix drops the stray argument. These
# tests pin the behavior that must survive the fix; they pass before and
# after, guarding against a regression while the call is edited.

test_that("small-sample branch returns a valid p-value under parallel = 'no'", {
  set.seed(20260811)
  smalldat <- data.frame(
    blockF = factor(rep("b1", 12)),
    trtF = factor(rep(c(0L, 1L), each = 6)),
    Y1 = rnorm(12) + 1.5 * rep(c(0, 1), each = 6)
  )
  p <- pWilcox(
    fmla = Y1 ~ trtF | blockF, dat = smalldat,
    simthresh = 20, sims = 500, parallel = "no"
  )
  expect_true(is.numeric(p) && length(p) == 1)
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("the same data above simthresh takes the asymptotic branch and
          agrees in direction", {
  # Not a claim of numerical equality -- Monte Carlo permutation and the
  # Strasser-Weber asymptotic reference differ in shape by design. Both
  # must see the same strong effect. This documents the branch boundary
  # the paper's supplement must describe accurately (Fix 4 doc change).
  set.seed(20260811)
  smalldat <- data.frame(
    blockF = factor(rep("b1", 12)),
    trtF = factor(rep(c(0L, 1L), each = 6)),
    Y1 = rnorm(12) + 3 * rep(c(0, 1), each = 6)
  )
  p_mc <- pWilcox(
    fmla = Y1 ~ trtF | blockF, dat = smalldat,
    simthresh = 20, sims = 1000, parallel = "no"
  )
  p_asy <- pWilcox(
    fmla = Y1 ~ trtF | blockF, dat = smalldat,
    simthresh = 5, sims = 1000, parallel = "no"
  )
  expect_lt(p_mc, 0.05)
  expect_lt(p_asy, 0.05)
})
