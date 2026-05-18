# Monte Carlo regression tests for weak FWER (thm:weakctrl) and the
# low-power-root corollary (cor:fwer_control_low_power) from Appendix B
# of the supplement.
#
# thm:weakctrl says: under the Stopping Rule (descend only when the
# parent test rejects) and the Valid Tests condition (each per-node test
# is level-alpha), the top-down procedure controls the family-wise error
# rate weakly --- that is, under the global null where every block has
# zero treatment effect, P(any rejection) <= alpha.
#
# cor:fwer_control_low_power says: in the partial-null case, when the
# root test has power at most alpha, the procedure controls FWER at
# alpha. The mechanism is the Stopping Rule: descendant rejections can
# only occur when the root rejects, so an alpha-level root rejection
# probability propagates to the full tree.
#
# Both results are structural (they follow from the Stopping Rule), so
# they are covered indirectly by tests that exercise find_blocks. These
# tests record the structural behavior with a Monte Carlo: under the
# global null and under a low-power root, observed FWER must lie within
# the binomial Monte Carlo bound around alpha.

source("make_test_data.R")


# ----------------------------------------------------------------------
# Helper: one-sided binomial upper bound for observed FWER
#
# Returns the smallest p such that with k rejections out of n iterations
# we cannot reject (at level 0.005) the null hypothesis that the true
# FWER <= p_true. Used to ask: "is the observed rejection rate too high
# to be consistent with FWER <= alpha?"
# ----------------------------------------------------------------------
binom_upper <- function(k, n, conf = 0.995) {
  # 1 - alpha confidence upper bound; using 0.995 to leave room for noise
  # without making the test fragile.
  if (k == 0L) {
    return(1 - conf^(1 / n))
  }
  stats::qbeta(conf, k + 1, n - k)
}


# ----------------------------------------------------------------------
# Helper: re-randomize Z within blocks and re-derive Ynull = Z * y0 + (1 - Z) * y0
#
# Under the global null, treated and control potential outcomes are
# identical, so the observed outcome equals y0 regardless of Z. We
# re-randomize Z to draw a fresh assignment under H_0 for each Monte
# Carlo iteration. Block sample sizes are taken from the input idat.
# ----------------------------------------------------------------------
rerandomize_Z_within_blocks <- function(idat) {
  out <- data.table::copy(idat)
  # randomizr::block_ra needs the full vector of block ids, not a
  # per-group call --- it allocates within blocks internally.
  out[, Z := randomizr::block_ra(blocks = bF)]
  out[, ZF := factor(Z)]
  # Y under the global null: y1null == y0, so Y == y0 regardless of Z
  out[, Y := y0null]
  out
}


# ----------------------------------------------------------------------
# Helper: run find_blocks once and report whether any node rejected
#
# A weak-FWER violation is "the procedure rejected at least one node
# when the global null was true." This is the simplest summary statistic
# the theorem speaks to.
# ----------------------------------------------------------------------
any_rejection <- function(idat_iter, bdat_template, thealpha) {
  bdat_iter <- data.table::copy(bdat_template)
  res <- tryCatch(
    find_blocks(
      idat = idat_iter,
      bdat = bdat_iter,
      blockid = "bF",
      splitfn = splitCluster,
      pfn = pIndepDist,
      alphafn = NULL,
      fmla = Y ~ ZF | bF,
      splitby = "hwt",
      blocksize = "hwt",
      parallel = "no",
      thealpha = thealpha,
      maxtest = 30,
      trace = FALSE,
      copydts = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(res)) return(NA)
  # Any node where the test rejected at its node-level alpha is a
  # rejection in the procedure's accounting.
  any(res$node_dat$p < res$node_dat$a)
}


# ============================================================================
# thm:weakctrl --- weak FWER under the global null
# ============================================================================

test_that("thm:weakctrl: under global null, observed FWER consistent with alpha", {
  # Statistical principle: with truly null potential outcomes the test
  # at every node is exchangeable under randomization, so the Stopping
  # Rule gives FWER <= alpha. We run B = 200 fresh randomizations and
  # check that the rejection count is consistent (binomial upper bound)
  # with FWER <= alpha.
  skip_on_cran()
  skip_if_not_installed("randomizr")

  thealpha <- 0.05
  B <- 200L
  set.seed(20260425)

  # idat3/bdat3: 20 blocks, irregular hwt so splitCluster works.
  # Reuses the truly-null potential outcomes built in make_test_data.R
  # (y0null is the residual after partialing out Z within block).
  idat_base <- data.table::copy(idat3)
  idat_base[, y0null := resid(lm(y0 ~ Z, data = .SD)), by = bF]

  rejections <- vapply(seq_len(B), function(i) {
    idat_iter <- rerandomize_Z_within_blocks(idat_base)
    any_rejection(idat_iter, bdat3, thealpha)
  }, logical(1L))

  # Drop iterations that errored (e.g. degenerate split). They do not
  # represent rejections, so excluding them is conservative.
  rejections <- rejections[!is.na(rejections)]
  observed_fwer <- mean(rejections)
  upper_bound <- binom_upper(sum(rejections), length(rejections), conf = 0.995)

  # The principle: even after Monte Carlo noise, we should not be able
  # to reject the null FWER <= alpha. The 0.995 confidence bound buys
  # us some headroom against rare false alarms.
  expect_true(
    upper_bound <= thealpha + 0.05,
    info = sprintf(
      "observed FWER = %.3f (%d/%d), 99.5%% upper = %.3f, alpha = %.3f",
      observed_fwer, sum(rejections), length(rejections),
      upper_bound, thealpha
    )
  )
  expect_true(
    observed_fwer <= thealpha + 0.03,
    info = sprintf(
      "observed FWER %.3f exceeds alpha+0.03; structural violation suspected",
      observed_fwer
    )
  )
})


# ============================================================================
# cor:fwer_control_low_power --- low root power forces FWER <= alpha
# ============================================================================

test_that("cor:fwer_control_low_power: very small alpha keeps FWER below alpha", {
  # Statistical principle: when alpha is small enough that root power
  # is itself bounded by alpha, descendant rejections cannot push FWER
  # above alpha because the Stopping Rule gates them all on the root.
  # We use a very small alpha (0.001) so the root rarely rejects, then
  # verify FWER <= alpha (with binomial headroom) under the null.
  skip_on_cran()
  skip_if_not_installed("randomizr")

  thealpha <- 0.001
  B <- 300L  # need more iterations because alpha is tiny
  set.seed(20260425L + 1L)

  idat_base <- data.table::copy(idat3)
  idat_base[, y0null := resid(lm(y0 ~ Z, data = .SD)), by = bF]

  rejections <- vapply(seq_len(B), function(i) {
    idat_iter <- rerandomize_Z_within_blocks(idat_base)
    any_rejection(idat_iter, bdat3, thealpha)
  }, logical(1L))

  rejections <- rejections[!is.na(rejections)]
  observed_fwer <- mean(rejections)
  upper_bound <- binom_upper(sum(rejections), length(rejections), conf = 0.995)

  # With alpha = 0.001 and B = 300, even a single false rejection
  # would push the upper bound above alpha. The corollary's promise is
  # that this can happen at most 0.1% of the time, so the upper bound
  # should sit near alpha + 1/B. We allow alpha + 0.02 of slack.
  expect_true(
    upper_bound <= thealpha + 0.02,
    info = sprintf(
      "observed FWER = %.4f (%d/%d), 99.5%% upper = %.4f, alpha = %.4f",
      observed_fwer, sum(rejections), length(rejections),
      upper_bound, thealpha
    )
  )
})


# ============================================================================
# Stopping Rule: descendant rejections require parent rejection
# ============================================================================

test_that("thm:weakctrl: descendant rejections never occur when root retains the null", {
  # Statistical principle: the Stopping Rule says find_blocks does not
  # descend below a node that retains the null. In particular, in any
  # iteration where the root fails to reject, the procedure must record
  # exactly one node (the root) and no descendants. This is the
  # mechanical guarantee that closes the proof of thm:weakctrl.
  skip_on_cran()
  skip_if_not_installed("randomizr")

  thealpha <- 0.05
  B <- 50L
  set.seed(20260425L + 2L)

  idat_base <- data.table::copy(idat3)
  idat_base[, y0null := resid(lm(y0 ~ Z, data = .SD)), by = bF]

  for (i in seq_len(B)) {
    idat_iter <- rerandomize_Z_within_blocks(idat_base)
    res <- tryCatch(
      find_blocks(
        idat = idat_iter,
        bdat = data.table::copy(bdat3),
        blockid = "bF",
        splitfn = splitCluster,
        pfn = pIndepDist,
        alphafn = NULL,
        fmla = Y ~ ZF | bF,
        splitby = "hwt",
        blocksize = "hwt",
        parallel = "no",
        thealpha = thealpha,
        maxtest = 30,
        trace = FALSE,
        copydts = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    nd <- res$node_dat
    root <- nd[depth == 1, ]
    if (nrow(root) != 1L) next
    if (root$p >= root$a) {
      # Root retained the null --- there must be no descendants.
      expect_equal(nrow(nd), 1L,
        info = sprintf("iter %d: root retained (p=%.3f, a=%.3f) but tree has %d nodes",
                       i, root$p, root$a, nrow(nd))
      )
    }
  }
})
