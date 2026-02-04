# Tests for the alphafn interface in find_blocks()
#
# Purpose: Verify that adding a `depth` parameter to the alphafn interface
# does not change the behavior of find_blocks() with the existing alpha
# functions (alpha_investing, alpha_saffron, alpha_addis).
#
# These tests record the current behavior and then check that it is
# preserved after the interface change.

source("make_test_data.R")

# --- Helper: run find_blocks with an alpha function on the test data ---
# Uses idat3/bdat3 (20 blocks with varying sizes/hwt) from make_test_data.R.
# The 10-block idat/bdat has constant hwt, which triggers the
# stop_splitby_constant guard in find_blocks.
run_fb_with_alpha <- function(alphafn, seed = 54321) {
  set.seed(seed)
  find_blocks(
    idat = idat3,
    bdat = data.table::copy(bdat3),
    blockid = "bF",
    splitfn = splitCluster,
    pfn = pOneway,
    alphafn = alphafn,
    fmla = Y ~ ZF | bF,
    splitby = "hwt",
    parallel = "no",
    thealpha = 0.05,
    thew0 = 0.049,
    maxtest = 15,
    trace = FALSE
  )
}

# --- Structural properties that must survive the interface change ---

test_that("find_blocks with alpha_investing returns expected structure", {
  res <- run_fb_with_alpha(alpha_investing)

  # Must return a list with bdat and node_dat
  expect_true(is.list(res))
  expect_true("bdat" %in% names(res))
  expect_true("node_dat" %in% names(res))

  # node_dat must have the columns find_blocks creates
  expect_true(all(c("parent", "p", "a", "depth", "nodesize", "nodenum") %in% names(res$node_dat)))

  # Root node (depth 1) always gets nominal alpha
  root <- res$node_dat[depth == 1]
  expect_equal(nrow(root), 1)
  expect_equal(root$a, 0.05)

  # Alpha values must all be positive
  # (alpha_investing can produce alphas above nominal when wealth accumulates,
  # so we only check positivity here)
  expect_true(all(res$node_dat$a > 0))

  # Every node must have a valid parent (0 for root, positive integer otherwise)
  expect_true(all(res$node_dat$parent >= 0))
  expect_equal(sum(res$node_dat$parent == 0), 1) # exactly one root

  # Depths must be sequential integers starting at 1
  expect_equal(min(res$node_dat$depth), 1)
  expect_true(all(diff(sort(unique(res$node_dat$depth))) == 1))
})

test_that("find_blocks with alpha_saffron returns expected structure", {
  res <- run_fb_with_alpha(alpha_saffron)

  expect_true(is.list(res))
  expect_true("bdat" %in% names(res))
  expect_true("node_dat" %in% names(res))

  root <- res$node_dat[depth == 1]
  expect_equal(nrow(root), 1)
  expect_equal(root$a, 0.05)

  expect_true(all(res$node_dat$a > 0))
  expect_equal(sum(res$node_dat$parent == 0), 1)
})

test_that("find_blocks with alpha_addis returns expected structure", {
  res <- run_fb_with_alpha(alpha_addis)

  expect_true(is.list(res))
  expect_true("bdat" %in% names(res))
  expect_true("node_dat" %in% names(res))

  root <- res$node_dat[depth == 1]
  expect_equal(nrow(root), 1)
  expect_equal(root$a, 0.05)

  expect_true(all(res$node_dat$a > 0))
  expect_equal(sum(res$node_dat$parent == 0), 1)
})

test_that("find_blocks with alphafn=NULL still works (baseline)", {
  # The no-alpha-adjustment case should also survive the change.
  # Uses idat3/bdat3 for consistency with the other tests.
  set.seed(54321)
  res <- find_blocks(
    idat = idat3,
    bdat = data.table::copy(bdat3),
    blockid = "bF",
    splitfn = splitCluster,
    pfn = pOneway,
    alphafn = NULL,
    fmla = Y ~ ZF | bF,
    splitby = "hwt",
    parallel = "no",
    thealpha = 0.05,
    maxtest = 15,
    trace = FALSE
  )

  expect_true(is.list(res))
  # When alphafn is NULL, all alpha values should be exactly thealpha
  expect_true(all(res$node_dat$a == 0.05))
})

test_that("existing alpha functions accept and ignore a depth argument", {
  # After the interface change, find_blocks will pass depth= to alphafn.
  # The existing functions must tolerate this extra argument.
  # This test will FAIL before we update the existing alpha functions,
  # which is the point — it defines what needs to change.

  pvals <- c(0.01, 0.04, 0.12)
  batches <- c(1, 2, 2)
  nodesizes <- c(100, 50, 50)
  depths <- c(1L, 2L, 2L)

  # First record results WITHOUT depth:
  res_inv_no_depth <- alpha_investing(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049)
  res_saf_no_depth <- alpha_saffron(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049)
  res_add_no_depth <- alpha_addis(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049)

  # These should be numeric vectors of the same length as input
  expect_length(res_inv_no_depth, 3)
  expect_length(res_saf_no_depth, 3)
  expect_length(res_add_no_depth, 3)

  # After we add `depth` (or `...`) to the function signatures, calling
  # with depth= should produce identical results because depth is unused
  # by these functions. This block will error until the signatures are updated:
  expect_no_error({
    res_inv_with_depth <- alpha_investing(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049, depth = depths)
  })
  expect_no_error({
    res_saf_with_depth <- alpha_saffron(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049, depth = depths)
  })
  expect_no_error({
    res_add_with_depth <- alpha_addis(pvals, batches, nodesizes, thealpha = 0.05, thew0 = 0.049, depth = depths)
  })

  # Results must be identical with or without depth
  expect_equal(res_inv_with_depth, res_inv_no_depth)
  expect_equal(res_saf_with_depth, res_saf_no_depth)
  expect_equal(res_add_with_depth, res_add_no_depth)
})
