# Tests for budget-weighted adaptive alpha (Theorem B.3 in supplement)
#
# The budget framework allocates a fraction w_ell of the total alpha
# budget to each depth. The adjusted test size is:
#   alpha_ell = min(alpha, w_ell * alpha / G_ell)
# where G_ell is the error load at depth ell. The constraint
# sum(w_ell) <= 1 guarantees FWER <= alpha by the union bound.
#
# Budget weights are needed when:
# - The tree is irregular (telescoping fails)
# - Branch pruning is used (denominators are data-dependent)
# - The user wants explicit control over the depth-wise allocation
#
# The default (budget_weights = NULL) preserves the current behavior:
# alpha_ell = alpha / G_ell, which is the tightest bound for regular
# k-ary trees via the telescoping argument.


# ============================================================================
# Helper: build a regular k-ary tree
# ============================================================================

make_regular_tree <- function(k, max_depth, N_total) {
  nodes <- data.frame(
    nodenum = 1L,
    parent = 0L,
    depth = 1L,
    nodesize = N_total
  )
  next_id <- 2L
  for (d in 2:max_depth) {
    parents_at_prev <- nodes$nodenum[nodes$depth == d - 1L]
    n_at_d <- N_total / k^(d - 1L)
    for (p in parents_at_prev) {
      for (j in seq_len(k)) {
        nodes <- rbind(nodes, data.frame(
          nodenum = next_id,
          parent = p,
          depth = d,
          nodesize = n_at_d
        ))
        next_id <- next_id + 1L
      }
    }
  }
  nodes
}


# ============================================================================
# Backward compatibility: NULL gives identical behavior to current code
# ============================================================================

test_that("budget_weights = NULL preserves telescoping behavior", {
  # A 3-ary tree with high power (error load > 1) where adjustment matters
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  # Current behavior (no budget_weights argument)
  alphas_old <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )

  # New behavior with budget_weights = NULL should be identical
  alphas_new <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05,
    budget_weights = NULL
  )

  expect_equal(as.numeric(alphas_new), as.numeric(alphas_old),
    info = "NULL budget_weights must reproduce the telescoping formula exactly"
  )
})


# ============================================================================
# Equal weights: w_ell = 1/(L-1) for all depths >= 2
# ============================================================================

test_that("equal budget weights split budget evenly across depths", {
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)

  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05,
    budget_weights = "equal"
  )

  # With L=3, there are 2 non-root depths, so w_ell = 1/2 each
  # alpha_ell = min(0.05, 0.5 * 0.05 / G_ell)
  # The key statistical property: sum of (w_ell * alpha) across depths = alpha
  # so FWER <= alpha by union bound regardless of null configuration
  el <- compute_error_load(node_dat = nd, delta_hat = 0.5)

  for (d in 2:3) {
    expected <- min(0.05, 0.5 * 0.05 / el$G[d])
    expect_equal(unname(alphas[d]), expected, tolerance = 1e-3,
      info = paste("Depth", d, "should use w = 1/2")
    )
  }
})


# ============================================================================
# Proportional weights: w_ell = G_ell / sum(G)
# ============================================================================

test_that("proportional weights give uniform adjusted alpha across depths", {
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)

  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05,
    budget_weights = "proportional"
  )

  # With proportional weights, alpha_ell = (G_ell / sum(G)) * alpha / G_ell
  # = alpha / sum(G) for all depths. All non-root depths get the same alpha.
  # This is a global Bonferroni correction by the total error load.
  # The key statistical property: uniform alpha across depths.
  expect_equal(unname(alphas[2]), unname(alphas[3]), tolerance = 1e-10,
    info = "Proportional weights should give the same alpha at every depth"
  )
  # And the common value should be less than nominal
  expect_true(alphas[2] < 0.05,
    info = "Proportional alpha should be reduced when error load > 1"
  )
})


# ============================================================================
# Custom numeric weights
# ============================================================================

test_that("custom numeric weights are applied correctly", {
  nd <- make_regular_tree(k = 3, max_depth = 4, N_total = 3000)

  # Front-load: give 60% to depth 2, 30% to depth 3, 10% to depth 4
  w <- c(0.6, 0.3, 0.1)

  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05,
    budget_weights = w
  )

  # The key statistical property: higher weight -> higher (less conservative)
  # alpha at that depth. Depth 2 (w=0.6) should have higher alpha than
  # depth 3 (w=0.3), which should have higher alpha than depth 4 (w=0.1).
  expect_true(alphas[2] > alphas[3],
    info = "Higher weight at depth 2 should give higher alpha"
  )
  expect_true(alphas[3] > alphas[4],
    info = "Higher weight at depth 3 should give higher alpha than depth 4"
  )
  # All should be below nominal
  expect_true(all(alphas[2:4] <= 0.05),
    info = "All adjusted alphas must be <= nominal"
  )
})

test_that("custom weights that sum to > 1 are rejected", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  # These weights sum to 1.2 > 1, violating the FWER guarantee
  expect_error(
    compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, budget_weights = c(0.7, 0.5)
    ),
    "sum to at most 1",
    info = "Weights summing to > 1 break the FWER guarantee"
  )
})

test_that("negative weights are rejected", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  expect_error(
    compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, budget_weights = c(0.8, -0.1)
    ),
    "non-negative",
    info = "Negative weights have no statistical meaning"
  )
})

test_that("wrong-length weight vector is rejected", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  # L = 3, so we need L-1 = 2 weights, not 3
  expect_error(
    compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, budget_weights = c(0.3, 0.3, 0.3)
    ),
    "length",
    info = "Weight vector must have one entry per non-root depth"
  )
})


# ============================================================================
# FWER guarantee: all adjusted alphas are <= nominal alpha
# ============================================================================

test_that("adjusted alphas never exceed nominal alpha", {
  nd <- make_regular_tree(k = 5, max_depth = 3, N_total = 5000)

  for (bw in list(NULL, "equal", "proportional", c(0.7, 0.3))) {
    alphas <- compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, thealpha = 0.05,
      budget_weights = bw
    )
    expect_true(all(alphas <= 0.05 + 1e-10),
      info = paste("All alphas must be <= nominal for budget_weights =",
                   deparse(bw))
    )
  }
})


# ============================================================================
# Budget weights are more conservative than telescoping for regular trees
# ============================================================================

test_that("budget weights are more conservative than telescoping on regular trees", {
  # On a regular tree, the telescoping formula (NULL) gives the tightest
  # bound. Any budget-weighted scheme must be at least as conservative.
  # This is the statistical content: the structural identity
  # m_ell + e_ell = k * m_{ell-1} lets telescoping avoid the union-bound
  # tax that budget weights must pay.
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)

  alphas_telescoping <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, budget_weights = NULL
  )

  alphas_equal <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, budget_weights = "equal"
  )

  # Budget-weighted alphas should be <= telescoping alphas at every depth
  # (more conservative = smaller alpha = fewer rejections)
  expect_true(all(alphas_equal[2:3] <= alphas_telescoping[2:3] + 1e-10),
    info = "Budget weights cannot beat telescoping on regular trees"
  )
})


# ============================================================================
# Natural gating regime: budget weights are irrelevant
# ============================================================================

test_that("budget weights have no effect when natural gating suffices", {
  # Small effect size so error load < 1 (natural gating regime)
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 300)

  alphas_null <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.1, budget_weights = NULL
  )

  alphas_equal <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.1, budget_weights = "equal"
  )

  # When error load < 1, both should return nominal alpha everywhere
  # because natural gating provides sufficient protection
  expect_equal(as.numeric(alphas_null), as.numeric(alphas_equal),
    info = "In natural gating regime, budget weights should not matter"
  )
})


# ============================================================================
# Factory function propagation: alpha_adaptive_tree passes budget_weights
# ============================================================================

test_that("alpha_adaptive_tree passes budget_weights through", {
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)

  # Create factory with equal weights
  fn <- alpha_adaptive_tree(
    node_dat = nd, delta_hat = 0.5, budget_weights = "equal"
  )

  # Call the factory function with depth info
  result <- fn(pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
               thealpha = 0.05, depth = c(1, 2, 3))

  # Compare with direct computation
  expected <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, budget_weights = "equal"
  )

  expect_equal(result, expected[c(1, 2, 3)],
    info = "Factory function should produce same alphas as direct computation"
  )
})
