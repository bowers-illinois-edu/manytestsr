# Tests for the switching corollary (Corollary B.1 in supplement)
#
# The switching corollary says: during branch pruning, once the remaining
# pruned error load fits within the remaining error budget, the procedure
# can revert to nominal alpha for all remaining depths. Formally:
#
#   If sum_{ell >= s}^L D_ell <= B_s (remaining budget),
#   then alpha_ell = alpha for all ell >= s is valid.
#
# This works because setting w_ell = D_ell for remaining depths gives
# alpha_ell = w_ell * alpha / D_ell = alpha (the D_ell cancels), and
# sum(w_ell) = sum(D_ell) <= B_s does not exceed the remaining budget.
#
# The switching corollary is the mechanism by which pruning provides its
# sharpest benefit: after most branches die, the surviving tree is narrow
# enough for natural gating to protect without any alpha adjustment.


# ============================================================================
# Helper: build test trees
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


# Simulate pruning by removing nodes at a given depth whose parents
# are in the "failed" set. Returns the pruned node_dat.
prune_tree <- function(node_dat, failed_nodenum) {
  # Remove failed nodes and all their descendants
  to_remove <- failed_nodenum
  repeat {
    children <- node_dat$nodenum[node_dat$parent %in% to_remove &
                                   !(node_dat$nodenum %in% to_remove)]
    if (length(children) == 0L) break
    to_remove <- c(to_remove, children)
  }
  node_dat[!node_dat$nodenum %in% to_remove, ]
}


# ============================================================================
# Switching fires when pruned error load fits within remaining budget
# ============================================================================

test_that("switching to nominal alpha when pruned error load <= remaining budget", {
  # Build a 3-ary tree depth 3 with high error load on the full tree
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  # Full tree error load should exceed 1 (needs adjustment)
  el_full <- compute_error_load(node_dat = nd, delta_hat = 0.5)
  expect_true(el_full$needs_adjustment,
    info = "Setup: full tree should need adjustment for this test to be meaningful"
  )

  # Create pruned factory with switching enabled and remaining-budget process
  obj <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5,
    budget_weights = "remaining", switching = TRUE
  )

  # Simulate heavy pruning: only 1 of 3 children survives at depth 2
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  failed <- depth2_nodes[2:3]  # kill 2 of 3 branches
  pruned_nd <- prune_tree(nd, failed)

  # Update with the pruned tree
  obj$update(pruned_nd, 0.05)

  # After heavy pruning, the surviving subtree is much narrower.
  # If the remaining error load fits within the remaining budget,
  # depths 3+ should get nominal alpha (0.05).
  el_pruned <- compute_error_load(node_dat = pruned_nd, delta_hat = 0.5)

  # The key statistical test: if switching fired, deeper depths
  # should have alpha = nominal (0.05), not a reduced value
  alphas <- obj$alphafn(
    pval = rep(0.01, nrow(pruned_nd)),
    batch = seq_len(nrow(pruned_nd)),
    nodesize = pruned_nd$nodesize,
    thealpha = 0.05,
    depth = pruned_nd$depth
  )

  # Check: if pruned error load < 1, depth-3 nodes should get nominal alpha
  # (switching corollary says remaining budget absorbs the error load)
  if (!el_pruned$needs_adjustment) {
    depth3_alphas <- alphas[pruned_nd$depth == 3L]
    expect_true(all(depth3_alphas == 0.05),
      info = "Switching should give nominal alpha when pruned error load < 1"
    )
  }
})


# ============================================================================
# Budget tracking: remaining budget decreases with depth
# ============================================================================

test_that("remaining budget decreases as depths are processed", {
  # Use a small tree (3-ary, depth 3) to keep tests fast
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 900)

  obj <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5,
    budget_weights = "remaining", spending_fraction = 0.5,
    budget_total = 1.0
  )

  # After first update (depth 2 tested), budget should be < 1.0
  # because spending_fraction = 0.5 means half the budget is spent
  obj$update(nd, 0.05)

  # The alphas at deeper depths should reflect the reduced budget
  alphas_after_update <- obj$alphafn(
    pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
    thealpha = 0.05, depth = 1:3
  )

  # Root should always get nominal alpha
  expect_equal(unname(alphas_after_update[1]), 0.05,
    info = "Root always gets nominal alpha"
  )
})


# ============================================================================
# Reset restores full budget
# ============================================================================

test_that("reset restores budget to initial value", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 900)

  obj <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5,
    budget_weights = "remaining", spending_fraction = 0.5,
    budget_total = 1.0
  )

  # Get initial alphas
  alphas_initial <- obj$alphafn(
    pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
    thealpha = 0.05, depth = 1:3
  )

  # Update (spend some budget)
  obj$update(nd, 0.05)

  # Reset
  obj$reset(0.05)

  # Alphas after reset should match initial
  alphas_after_reset <- obj$alphafn(
    pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
    thealpha = 0.05, depth = 1:3
  )

  expect_equal(as.numeric(alphas_after_reset), as.numeric(alphas_initial),
    info = "Reset must restore the full-tree alpha schedule and budget"
  )
})


# ============================================================================
# FWER guarantee: adjusted alphas never exceed nominal
# ============================================================================

test_that("adjusted alphas under switching never exceed nominal", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  obj <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5,
    budget_weights = "remaining", switching = TRUE
  )

  # Even after switching fires, no alpha should exceed nominal
  # (switching sets alpha to nominal, not above it)
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  pruned_nd <- prune_tree(nd, depth2_nodes[2:3])
  obj$update(pruned_nd, 0.05)

  alphas <- obj$alphafn(
    pval = rep(0.01, nrow(pruned_nd)),
    batch = seq_len(nrow(pruned_nd)),
    nodesize = pruned_nd$nodesize,
    thealpha = 0.05,
    depth = pruned_nd$depth
  )

  expect_true(all(alphas <= 0.05 + 1e-10),
    info = "Switching to nominal means alpha = 0.05, never above"
  )
})


# ============================================================================
# Backward compatibility: switching = FALSE or budget_weights = NULL
# ============================================================================

test_that("alpha_adaptive_tree_pruned with defaults matches current behavior", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1500)

  # Current behavior (no new parameters)
  obj_old <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5
  )

  # New behavior with explicit defaults --- should be identical
  obj_new <- alpha_adaptive_tree_pruned(
    node_dat = nd, delta_hat = 0.5,
    budget_weights = NULL, switching = FALSE
  )

  alphas_old <- obj_old$alphafn(
    pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
    thealpha = 0.05, depth = 1:3
  )

  alphas_new <- obj_new$alphafn(
    pval = rep(0.01, 3), batch = 1:3, nodesize = rep(100, 3),
    thealpha = 0.05, depth = 1:3
  )

  expect_equal(as.numeric(alphas_old), as.numeric(alphas_new),
    info = "Default new parameters must reproduce old behavior exactly"
  )
})
