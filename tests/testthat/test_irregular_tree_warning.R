# Tests for the irregular-tree FWER guard rail.
#
# The per-depth telescoping formula alpha_ell = alpha / G_ell controls
# FWER exactly for *regular* k-ary trees (Theorem B.2), but it can
# under-protect on irregular trees with sum_G > 1 because the path
# powers in different branches can be unequal in ways that the
# depth-wise sum hides. The supplement's irregular-tree memo
# (Theory/irregular_tree_branch_pruning_analysis.md) gives a
# counterexample where this happens.
#
# The supported safe path on irregular trees is the budget-weighted
# formulation (Theorem B.3): pass budget_weights = "equal" or
# "proportional" or a custom numeric vector. The warning here exists
# so a user accidentally relying on telescoping in the unsafe regime
# notices --- the math is fine when sum_G <= 1 (natural gating) so we
# stay quiet there.
#
# Trigger conditions: tree is irregular AND sum_G > 1 AND budget_weights
# is NULL. All three must hold.

# ============================================================================
# Helpers --- regular and irregular trees
# ============================================================================

.regular_tree <- function(k, max_depth, N_total) {
  nd <- data.frame(
    nodenum = 1L, parent = 0L, depth = 1L, nodesize = N_total
  )
  next_id <- 2L
  for (d in 2:max_depth) {
    parents_at_prev <- nd$nodenum[nd$depth == d - 1L]
    n_at_d <- N_total / k^(d - 1L)
    for (p in parents_at_prev) {
      for (j in seq_len(k)) {
        nd <- rbind(nd, data.frame(nodenum = next_id, parent = p,
                                   depth = d, nodesize = n_at_d))
        next_id <- next_id + 1L
      }
    }
  }
  nd
}

.irregular_tree_high_load <- function() {
  # Skewed tree: root --> 2 children at depth 2 --> {3, 1} children at
  # depth 3. Same shape as the irregular-tree counterexample. We give
  # it a large root so theta is high enough to push sum_G above 1.
  data.frame(
    nodenum  = 1:7,
    parent   = c(0L, 1L, 1L, 2L, 2L, 2L, 3L),
    depth    = c(1L, 2L, 2L, 3L, 3L, 3L, 3L),
    nodesize = c(2000, 1000, 1000, 333, 333, 333, 1000)
  )
}

.irregular_tree_low_load <- function() {
  # Same shape as the high-load tree but with a small root, so power
  # is low everywhere and natural gating (sum_G <= 1) holds.
  data.frame(
    nodenum  = 1:7,
    parent   = c(0L, 1L, 1L, 2L, 2L, 2L, 3L),
    depth    = c(1L, 2L, 2L, 3L, 3L, 3L, 3L),
    nodesize = c(10, 5, 5, 2, 2, 2, 5)
  )
}


# ============================================================================
# Warning fires only on the (irregular & sum_G > 1 & telescoping) path
# ============================================================================

test_that("warning when irregular tree has high error load and budget_weights = NULL", {
  # Statistical principle: this is the regime where telescoping under-
  # protects (per the irregular-tree memo). The user should be told.
  nd <- .irregular_tree_high_load()
  el <- compute_error_load(node_dat = nd, delta_hat = 0.5)
  expect_true(el$needs_adjustment,
    info = "Setup: tree must have sum_G > 1 for the warning regime"
  )
  expect_warning(
    compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5),
    "irregular"
  )
})

test_that("no warning when irregular tree has low error load (natural gating)", {
  # Statistical principle: when sum_G <= 1, every adjustment formula
  # collapses to nominal alpha (Theorem B.1). The telescoping question
  # never arises, so the warning is suppressed.
  nd <- .irregular_tree_low_load()
  el <- compute_error_load(node_dat = nd, delta_hat = 0.5)
  expect_false(el$needs_adjustment,
    info = "Setup: tree must have sum_G <= 1 for natural gating"
  )
  expect_no_warning(
    compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5)
  )
})

test_that("no warning when irregular tree uses budget_weights", {
  # Statistical principle: budget_weights = 'equal' / 'proportional' /
  # numeric vector is the supported irregular-tree path (Theorem B.3),
  # so we stay quiet when the user has selected one.
  nd <- .irregular_tree_high_load()
  expect_no_warning(
    compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5,
                                 budget_weights = "equal")
  )
  expect_no_warning(
    compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5,
                                 budget_weights = "proportional")
  )
})

test_that("no warning when tree is regular (telescoping is exact)", {
  # Statistical principle: telescoping is *the* tight bound for regular
  # k-ary trees (Theorem B.2). No warning needed.
  nd <- .regular_tree(k = 3, max_depth = 3, N_total = 1500)
  el <- compute_error_load(node_dat = nd, delta_hat = 0.5)
  expect_true(el$needs_adjustment,
    info = "Setup: regular tree should have sum_G > 1 here too"
  )
  expect_no_warning(
    compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5)
  )
})

test_that("alpha_adaptive_tree factory propagates the warning", {
  # The factory wraps compute_adaptive_alphas_tree, so the warning
  # should fire on the first cache miss --- typically at the first call
  # to the closure.
  nd <- .irregular_tree_high_load()
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)
  expect_warning(
    fn(pval = 0.01, batch = 1, nodesize = 100,
       thealpha = 0.05, thew0 = 0.049, depth = 2L),
    "irregular"
  )
})
