# Tests for withdrawing the pruned-schedule guarantee (referee-correction
# Fix 2, see FIX_PLAN.md). Written 2026-08-11 BEFORE the implementation fix.
#
# The statistical principle: the pruned-load schedule
# alpha_ell = min(alpha, w_ell * alpha / D_ell) does NOT control the FWER.
# Exact counterexample (no Monte Carlo error): a 3-depth tree with a weak
# root, 30 moderate-power parents, and 20 children per parent reaches
# FWER = 0.063 at alpha = 0.05 with valid tests, conservative power, and
# the intersection structure all holding. The mechanism: D_ell charges each
# exposed null its estimated path power (an unconditional reach
# probability), but conditional on its parent's rejection the null is a
# full test; when path powers are small the count of exposed nulls exceeds
# the charged load and the schedule over-spends. The switching rule
# inherits the same failure. Until a replacement schedule with a proof
# ships, the constructor must warn so no user mistakes this for a
# guarantee -- while still computing, so existing pipelines (and the
# submitted paper's numbers, tag pre-referee-fixes) stay reproducible.

# Small regular tree with headcount node sizes: root 900, three children
# of 300, nine grandchildren of 100.
make_k3_tree <- function() {
  nodes <- data.frame(nodenum = 1L, parent = 0L, depth = 1L, nodesize = 900)
  id <- 2L
  for (p in 1:1) {
    for (j in 1:3) {
      nodes <- rbind(nodes, data.frame(
        nodenum = id, parent = 1L, depth = 2L, nodesize = 300
      ))
      id <- id + 1L
    }
  }
  for (p in 2:4) {
    for (j in 1:3) {
      nodes <- rbind(nodes, data.frame(
        nodenum = id, parent = p, depth = 3L, nodesize = 100
      ))
      id <- id + 1L
    }
  }
  nodes
}

test_that("the pruned-load schedule warns that its guarantee does not hold", {
  tree <- make_k3_tree()
  expect_warning(
    alpha_adaptive_tree_pruned(node_dat = tree, delta_hat = 0.2),
    regexp = "does not hold"
  )
})

test_that("the switching path carries the same warning", {
  # cor:switching re-spends to nominal alpha exactly when pruning has made
  # the realized load small -- which the counterexample shows is when the
  # exposed-null count most exceeds the charged load. It inherits the
  # falsified argument and must not be presented as safe.
  tree <- make_k3_tree()
  expect_warning(
    alpha_adaptive_tree_pruned(
      node_dat = tree, delta_hat = 0.2,
      budget_weights = "depth-sequential", switching = TRUE
    ),
    regexp = "does not hold"
  )
})

test_that("the pruned constructor still returns a working interface", {
  # Withdrawal, not removal: the warning must not break replication of the
  # submitted results. The three-part interface survives.
  tree <- make_k3_tree()
  obj <- suppressWarnings(
    alpha_adaptive_tree_pruned(
      node_dat = tree, delta_hat = 0.2,
      budget_weights = "depth-sequential", switching = TRUE
    )
  )
  expect_true(is.list(obj))
  expect_true(all(c("alphafn", "update", "reset") %in% names(obj)))
  expect_true(is.function(obj$alphafn))
})

test_that("the static budget schedule (Theorem 4) does not warn", {
  # The sound alternative must stay warning-free so the message
  # differentiates rather than blankets: static weights over the full-tree
  # load survived the same 45,000-configuration exact search that
  # falsified the pruned schedule.
  tree <- make_k3_tree()
  expect_silent(
    compute_adaptive_alphas_tree(
      node_dat = tree, delta_hat = 0.2, thealpha = 0.05,
      budget_weights = "equal"
    )
  )
})
