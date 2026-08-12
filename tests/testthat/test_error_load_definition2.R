# Tests for the Definition-2 error load (referee-correction Fix 1, see
# FIX_PLAN.md). Written 2026-08-11 BEFORE the implementation fix: each test
# encodes the paper's Definition 2 and the two-tailed power model, and FAILS
# against the pre-fix code, which (a) uses one-tailed power, (b) multiplies
# each node's own theta back into the load, (c) sums from depth 1, and
# (d) silently accepts weight-scale node sizes.
#
# The statistical principle: the error load G_ell is the expected number of
# depth-ell nodes the procedure REACHES (sum over depth-ell nodes of the
# product of ancestor rejection probabilities), not the expected number it
# rejects. The FWER bound multiplies each reached null by its test level, so
# including the node's own theta double-counts the rejection step, and the
# root (depth 1, always tested) contributes no reach term. Understating the
# load flips needs_adjustment to FALSE when the truth is "adjust" -- the
# anti-conservative direction the AOAS referee identified.

# Independent transcription of the corrected power model: two-tailed
# level-alpha z-test power. At delta_hat = 0 this equals the SIZE alpha
# (both tails), which is the referee's floor-doubling point. Kept separate
# from package internals so a shared bug cannot hide.
two_tailed_power <- function(delta_hat, n, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  pnorm(delta_hat * sqrt(n) - z) + pnorm(-delta_hat * sqrt(n) - z)
}

# Hand-computable fixture: root (400 units) with two children of 200; the
# second child has two children of 100. Node sizes are headcounts.
#   Definition 2: G_2 = 2 * theta_root; G_3 = 2 * theta_root * theta_3;
#   sum_G = 2 * theta_root * (1 + theta_3). No depth-1 term.
small_tree <- data.frame(
  nodenum = 1:5,
  parent = c(0L, 1L, 1L, 3L, 3L),
  depth = c(1L, 2L, 2L, 3L, 3L),
  nodesize = c(400, 200, 200, 100, 100)
)

test_that("power model has two tails: theta at delta_hat = 0 equals alpha", {
  el <- manytestsr:::.error_load_from_tree(
    small_tree,
    delta_hat = 0, z_crit = qnorm(1 - 0.05 / 2), thealpha = 0.05
  )
  # At zero effect the 'power' is the size of the two-sided test: 0.05.
  # The pre-fix one-tailed formula returns alpha/2 = 0.025 -- half the
  # answer exactly where the MDRC thetas sit (the referee's point 1d).
  expect_equal(el$node_detail$theta, rep(0.05, 5), tolerance = 1e-10)
})

test_that("sum_G matches Definition 2 on the hand-computed tree", {
  d <- 0.20
  el <- compute_error_load(node_dat = small_tree, delta_hat = d)
  t_root <- two_tailed_power(d, 400)
  t_3 <- two_tailed_power(d, 200)
  expected <- 2 * t_root + 2 * t_root * t_3
  expect_equal(el$sum_G, expected, tolerance = 1e-8)
})

test_that("a root-only tree carries zero error load", {
  # The root is tested at alpha whether or not we compute a load; the load
  # counts EXPOSURE below the root. Definition 2 sums from depth 2, so a
  # single-node tree has sum_G = 0. The pre-fix code returns theta_root.
  root_only <- data.frame(
    nodenum = 1L, parent = 0L, depth = 1L, nodesize = 500
  )
  el <- compute_error_load(node_dat = root_only, delta_hat = 0.3)
  expect_equal(el$sum_G, 0)
  expect_false(el$needs_adjustment)
})

test_that("the needs_adjustment gate agrees with the alpha denominators", {
  # A flat tree: root (300 units) with three children of 100. Choose
  # delta_hat so that the Definition-2 load 3 * theta_root just exceeds 1.
  # The pre-fix gate (own-theta formula) evaluates well below 1 here and
  # returns nominal alpha everywhere -- the internal disagreement the
  # referee flagged between R/alpha_adaptive.R:517 and :544: the quantity
  # deciding "no adjustment needed" must be the quantity that bounds the
  # FWER.
  flat_tree <- data.frame(
    nodenum = 1:4,
    parent = c(0L, 1L, 1L, 1L),
    depth = c(1L, 2L, 2L, 2L),
    nodesize = c(300, 100, 100, 100)
  )
  d <- 0.0985
  # Fixture validation, not the claim under test: Definition-2 load > 1.
  expect_gt(3 * two_tailed_power(d, 300), 1)
  alphas <- compute_adaptive_alphas_tree(
    node_dat = flat_tree, delta_hat = d, thealpha = 0.05
  )
  # With sum_G > 1 the schedule must tighten somewhere below the root.
  expect_true(any(alphas[-1] < 0.05))
})

test_that("weight-scale node sizes are refused, not silently used", {
  # find_blocks' default blocksize = "hwt" produces node sizes on the unit
  # interval (the supplement F.2 root printed 0.2254 for 1,268 people).
  # Power at a 'sample size' of 0.25 is meaningless; feeding it to
  # pnorm(delta*sqrt(n) - z) floors every theta near alpha/2 and produced
  # the paper's false all-clear. The power path must refuse counts below 1
  # with a message that names the trap.
  hwt_tree <- data.frame(
    nodenum = 1:3,
    parent = c(0L, 1L, 1L),
    depth = c(1L, 2L, 2L),
    nodesize = c(0.2254, 0.1030, 0.1224)
  )
  expect_error(
    compute_error_load(node_dat = hwt_tree, delta_hat = 0.2),
    regexp = "headcount|weight|nodesize"
  )
  expect_error(
    compute_adaptive_alphas_tree(
      node_dat = hwt_tree, delta_hat = 0.2, thealpha = 0.05
    ),
    regexp = "headcount|weight|nodesize"
  )
  # The pruned factory computes an error load at creation, so a
  # find_blocks(blocksize = "hwt") pipeline that feeds its node_dat here
  # fails loudly BEFORE any testing runs, not silently mid-analysis.
  expect_error(
    suppressWarnings(
      alpha_adaptive_tree_pruned(node_dat = hwt_tree, delta_hat = 0.2)
    ),
    regexp = "headcount|weight|nodesize"
  )
})
