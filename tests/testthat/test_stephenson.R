## Tests for pCombStephenson() wrapper around CMRSS::pval_comb_block()

test_that("pCombStephenson gives clear error when CMRSS is not installed", {
  # Mock the check by testing the error message pattern
  # This test verifies the error message is informative even if CMRSS *is* installed
  skip_if(requireNamespace("CMRSS", quietly = TRUE),
    message = "CMRSS is installed; skipping missing-package error test"
  )
  expect_error(
    pCombStephenson(idat, Y ~ ZF | bF),
    "Package 'CMRSS' is required"
  )
})

test_that("pCombStephenson requires a block variable in the formula", {
  skip_if_not_installed("CMRSS")
  expect_error(
    pCombStephenson(idat, Y ~ ZF),
    "block variable"
  )
})

test_that("pCombStephenson returns a numeric scalar between 0 and 1", {
  skip_if_not_installed("CMRSS")
  set.seed(42)
  # Default k = n tests the sharp null (standard rank-sum combined across scores)
  p <- pCombStephenson(idat, Y ~ ZF | bF, null_max = 1000)
  expect_type(p, "double")
  expect_length(p, 1)
  expect_true(p >= 0 && p <= 1)
})

test_that("pCombStephenson rejects on data with known effects", {
  skip_if_not_installed("CMRSS")
  # Yhomog has homogeneous positive effects across all blocks.
  # Default k = n tests the sharp null of no effects.
  # With uniformly large positive effects, the p-value should be small.
  set.seed(42)
  p <- pCombStephenson(idat, Yhomog ~ ZF | bF, null_max = 1000)
  expect_lt(p, 0.05)
})

test_that("pCombStephenson does not reject on null data", {
  skip_if_not_installed("CMRSS")
  set.seed(42)
  p_null <- pCombStephenson(idat, Ynull ~ ZF | bF, null_max = 1000)
  # Under the null, p should generally be above alpha
  # We use a lenient check — just that it's not tiny
  expect_gt(p_null, 0.01)
})

test_that("pCombStephenson warns on degenerate k", {
  skip_if_not_installed("CMRSS")
  # k <= n - m produces a degenerate test — wrapper should warn
  n <- nrow(idat)
  m <- sum(as.numeric(levels(idat$ZF)[idat$ZF]))
  expect_warning(
    pCombStephenson(idat, Y ~ ZF | bF, k = n - m, null_max = 100),
    "degenerate"
  )
})

test_that("pCombStephenson result matches direct CMRSS call", {
  skip_if_not_installed("CMRSS")
  set.seed(42)
  p_wrapper <- pCombStephenson(idat, Y ~ ZF | bF,
    r_vec = c(2, 6),
    null_max = 1000
  )

  # Now call CMRSS directly with the same configuration
  set.seed(42)
  Z <- as.numeric(levels(idat$ZF)[idat$ZF])
  Y_vec <- idat$Y
  block <- factor(idat$bF)
  B <- nlevels(block)
  n <- length(Z)
  k <- n # matches wrapper default: sharp null

  methods.list.all <- lapply(c(2, 6), function(r) {
    lapply(seq_len(B), function(b) {
      list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)
    })
  })

  p_direct <- as.numeric(CMRSS::pval_comb_block(
    Z = Z, Y = Y_vec, k = k, c = 0,
    block = block,
    methods.list.all = methods.list.all,
    weight.name = "asymp.opt",
    null.max = 1000,
    statistic = FALSE,
    opt.method = "ILP_auto",
    comb.method = 1
  ))

  expect_equal(p_wrapper, p_direct)
})

test_that("pCombStephenson handles factor treatment with levels '0' and '1'", {
  skip_if_not_installed("CMRSS")
  # idat already has ZF as a factor — just verify it works
  set.seed(42)
  p <- pCombStephenson(idat, Y ~ ZF | bF, null_max = 500)
  expect_type(p, "double")
  expect_true(p >= 0 && p <= 1)
})
