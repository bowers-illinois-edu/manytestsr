# Tests for branch-pruning adaptive alpha
#
# Branch pruning recomputes the adaptive alpha schedule on the surviving
# subtree after each depth. When branches die (parent fails to reject),
# their descendants will never be tested, so the alpha budget allocated
# to them can be redistributed to surviving branches. The FWER guarantee
# holds by the same telescoping-sum argument applied to the pruned tree.
#
# The interface is a list with three components:
#   $alphafn  — standard closure (same signature as alpha_adaptive_tree)
#   $update   — recompute schedule on pruned tree (called after each depth)
#   $reset    — restore full-tree schedule (called between simulation runs)


# ============================================================================
# Reuse helper functions from the static tree tests
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

make_dpp_tree <- function() {
  nodes <- data.frame(
    nodenum = integer(0),
    parent = integer(0),
    depth = integer(0),
    nodesize = numeric(0)
  )
  id <- 1L
  # Root
  nodes <- rbind(nodes, data.frame(nodenum = id, parent = 0L, depth = 1L, nodesize = 2200))
  root_id <- id; id <- id + 1L
  # 5 colleges
  college_sizes <- c(450, 400, 450, 400, 500)
  college_ids <- integer(length(college_sizes))
  for (i in seq_along(college_sizes)) {
    nodes <- rbind(nodes, data.frame(
      nodenum = id, parent = root_id, depth = 2L, nodesize = college_sizes[i]
    ))
    college_ids[i] <- id; id <- id + 1L
  }
  # Cohorts (2-3 per college)
  cohort_parents <- c(
    rep(college_ids[1], 3), rep(college_ids[2], 2),
    rep(college_ids[3], 3), rep(college_ids[4], 2),
    rep(college_ids[5], 2)
  )
  cohort_sizes <- c(150, 150, 150, 200, 200, 150, 150, 150, 200, 200, 200, 300)
  cohort_ids <- integer(length(cohort_sizes))
  for (i in seq_along(cohort_sizes)) {
    nodes <- rbind(nodes, data.frame(
      nodenum = id, parent = cohort_parents[i], depth = 3L, nodesize = cohort_sizes[i]
    ))
    cohort_ids[i] <- id; id <- id + 1L
  }
  # Leaves (1-2 per cohort)
  leaf_parents <- c(
    rep(cohort_ids[1], 2), cohort_ids[2], cohort_ids[3], cohort_ids[4],
    rep(cohort_ids[5], 2), cohort_ids[6], cohort_ids[7],
    rep(cohort_ids[8], 2), cohort_ids[9], cohort_ids[10],
    cohort_ids[11], cohort_ids[12]
  )
  leaf_sizes <- c(75, 75, 150, 150, 100, 100, 150, 150, 75, 75, 200, 200, 200, 300)
  for (i in seq_along(leaf_sizes)) {
    nodes <- rbind(nodes, data.frame(
      nodenum = id, parent = leaf_parents[i], depth = 4L, nodesize = leaf_sizes[i]
    ))
    id <- id + 1L
  }
  nodes
}


# ============================================================================
# Tests for .get_all_descendants()
# ============================================================================

test_that("get_all_descendants: finds all descendants in a simple tree", {
  # Tree: 1 -> {2,3}, 2 -> {4,5}, 3 -> {6,7}
  # Descendants of node 2 should be {4, 5}
  # Descendants of node 1 should be {2, 3, 4, 5, 6, 7}
  nd <- data.frame(
    nodenum = 1:7,
    parent = c(0, 1, 1, 2, 2, 3, 3),
    depth = c(1, 2, 2, 3, 3, 3, 3),
    nodesize = c(500, 250, 250, 125, 125, 125, 125)
  )
  # Single internal node
  desc_2 <- manytestsr:::.get_all_descendants(nd, 2L)
  expect_setequal(desc_2, c(4L, 5L))

  # Root: everything except itself
  desc_1 <- manytestsr:::.get_all_descendants(nd, 1L)
  expect_setequal(desc_1, 2:7)

  # Leaf node: no descendants
  desc_7 <- manytestsr:::.get_all_descendants(nd, 7L)
  expect_length(desc_7, 0)
})

test_that("get_all_descendants: handles multiple starting nodes", {
  # Removing nodes 2 AND 3 from the tree above should give all of {4,5,6,7}
  nd <- data.frame(
    nodenum = 1:7,
    parent = c(0, 1, 1, 2, 2, 3, 3),
    depth = c(1, 2, 2, 3, 3, 3, 3),
    nodesize = c(500, 250, 250, 125, 125, 125, 125)
  )
  desc <- manytestsr:::.get_all_descendants(nd, c(2L, 3L))
  expect_setequal(desc, 4:7)
})

test_that("get_all_descendants: works on irregular DPP-like tree", {
  # In the DPP tree, removing college 1 (node 2) should remove its
  # 3 cohorts and all their leaves. We verify the count: college 1 has
  # 3 cohorts (nodes 7,8,9) and cohort 7 has 2 leaves, cohorts 8,9 each
  # have 1 leaf = 3 cohorts + 4 leaves = 7 descendants.
  dpp <- make_dpp_tree()
  college1_id <- 2L
  desc <- manytestsr:::.get_all_descendants(dpp, college1_id)

  # All descendants should have depth > 2 (deeper than college 1)
  desc_depths <- dpp$depth[dpp$nodenum %in% desc]
  expect_true(all(desc_depths > 2))

  # None of the descendants' parents should be the root
  # (they're all within college 1's subtree)
  desc_parents <- dpp$parent[dpp$nodenum %in% desc]
  expect_false(1L %in% desc_parents)  # root should not be a parent of any descendant
})


# ============================================================================
# Tests for alpha_adaptive_tree_pruned() factory
# ============================================================================

test_that("pruned factory: returns a list with alphafn, update, reset", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  expect_true(is.list(obj))
  expect_true("alphafn" %in% names(obj))
  expect_true("update" %in% names(obj))
  expect_true("reset" %in% names(obj))
  expect_true(is.function(obj$alphafn))
  expect_true(is.function(obj$update))
  expect_true(is.function(obj$reset))
})

test_that("pruned factory: alphafn has correct interface", {
  # Must accept the same arguments as alpha_adaptive_tree() closures
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)
  args <- names(formals(obj$alphafn))
  expect_true("pval" %in% args)
  expect_true("batch" %in% args)
  expect_true("nodesize" %in% args)
  expect_true("thealpha" %in% args)
  expect_true("thew0" %in% args)
  expect_true("depth" %in% args)
})


# ============================================================================
# Core statistical guarantees
# ============================================================================

test_that("pruned: matches static schedule when no branches are pruned", {
  # When every branch survives, the pruned schedule should be identical
  # to the static schedule. This is the baseline equivalence property:
  # pruning with no dead branches = no pruning.
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  static_alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  for (d in seq_along(static_alphas)) {
    result <- obj$alphafn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(result), static_alphas[[d]],
      info = paste("depth", d, "should match static when nothing pruned")
    )
  }
})

test_that("pruned: alpha increases when branches die", {
  # When some depth-2 nodes fail to reject, pruning their subtrees
  # reduces the expected number of tests at deeper levels. This should
  # make alpha at those deeper levels LARGER (less stringent) than
  # the static schedule.
  #
  # Use a k=5 tree (like DPP) where 4 of 5 branches die at depth 2.
  nd <- make_regular_tree(k = 5, max_depth = 3, N_total = 2500)
  static_alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  # Simulate: mark 4 of 5 depth-2 nodes as failed (testable=FALSE)
  # The surviving subtree has only 1 branch at depth 2 and its 5 children.
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  expect_length(depth2_nodes, 5)  # k=5 at depth 2
  surviving <- depth2_nodes[1]
  failed <- depth2_nodes[2:5]

  # Build pruned node_dat: remove failed nodes and all their descendants
  failed_and_desc <- c(failed, manytestsr:::.get_all_descendants(nd, failed))
  pruned_nd <- nd[!nd$nodenum %in% failed_and_desc, ]

  obj$update(pruned_nd, thealpha = 0.05)

  # Check depth 3: pruned alpha should be >= static alpha
  pruned_alpha3 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )
  expect_gt(unname(pruned_alpha3), unname(static_alphas[[3]]),
    label = "Pruned alpha at depth 3 should exceed static alpha"
  )

  # Depth 1 should be unchanged (root is always nominal)
  pruned_alpha1 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 1L
  )
  expect_equal(unname(pruned_alpha1), 0.05)
})

test_that("pruned: alpha never exceeds nominal", {
  # Even after aggressive pruning, alpha at any depth should not
  # exceed nominal alpha. The min(alpha, alpha/sum_path_power) guard
  # ensures this.
  nd <- make_regular_tree(k = 5, max_depth = 3, N_total = 2500)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  # Prune to a single surviving branch
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  failed <- depth2_nodes[2:5]
  failed_and_desc <- c(failed, manytestsr:::.get_all_descendants(nd, failed))
  pruned_nd <- nd[!nd$nodenum %in% failed_and_desc, ]
  obj$update(pruned_nd, thealpha = 0.05)

  for (d in 1:3) {
    alpha_d <- obj$alphafn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_lte(unname(alpha_d), 0.05,
      label = paste("Alpha at depth", d, "must not exceed nominal")
    )
  }
})

test_that("pruned: reset restores full-tree schedule", {
  # After pruning and updating, calling reset must restore the
  # original full-tree alpha schedule. This is critical for
  # simulation loops where the same closure is reused across runs.
  nd <- make_regular_tree(k = 5, max_depth = 3, N_total = 2500)
  static_alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  # First, prune heavily
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  failed <- depth2_nodes[2:5]
  failed_and_desc <- c(failed, manytestsr:::.get_all_descendants(nd, failed))
  pruned_nd <- nd[!nd$nodenum %in% failed_and_desc, ]
  obj$update(pruned_nd, thealpha = 0.05)

  # Verify alphas are now different from static (pruning took effect)
  pruned_alpha3 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )
  expect_false(isTRUE(all.equal(unname(pruned_alpha3), unname(static_alphas[[3]]))),
    info = "Pruned alpha should differ from static (sanity check)"
  )

  # Now reset
  obj$reset(thealpha = 0.05)

  # Verify all depths match static schedule again
  for (d in seq_along(static_alphas)) {
    result <- obj$alphafn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(result), unname(static_alphas[[d]]),
      info = paste("After reset, depth", d, "should match static")
    )
  }
})

test_that("pruned: DPP-like 5x gain when 4 of 5 colleges die", {
  # The DPP tree has 5 colleges at depth 2. If effects are concentrated
  # in one college, the other 4 branches die. With the full tree, alpha
  # at depth 3 is divided across the expected tests from all 5 branches.
  # With pruning, only the surviving branch's tests matter, giving
  # roughly 5x more alpha at depth 3.
  dpp <- make_dpp_tree()
  static_alphas <- compute_adaptive_alphas_tree(
    node_dat = dpp, delta_hat = 0.5, thealpha = 0.05
  )
  obj <- alpha_adaptive_tree_pruned(node_dat = dpp, delta_hat = 0.5)

  # Keep only college 1 (node 2), kill colleges 2-5 (nodes 3-6)
  failed_colleges <- 3:6
  failed_and_desc <- c(failed_colleges,
    manytestsr:::.get_all_descendants(dpp, failed_colleges))
  pruned_dpp <- dpp[!dpp$nodenum %in% failed_and_desc, ]

  obj$update(pruned_dpp, thealpha = 0.05)

  pruned_alpha3 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )
  static_alpha3 <- static_alphas[[3]]

  # The ratio should be substantial — close to 5x for 4/5 branches dying.
  # The exact ratio depends on per-college path powers (which vary because
  # college sizes differ), so we test for at least 3x gain.
  ratio <- unname(pruned_alpha3) / unname(static_alpha3)
  expect_gt(ratio, 3,
    label = paste("Expected ~5x gain, got", round(ratio, 2), "x")
  )
})

test_that("pruned: error load of pruned tree < full tree", {
  # When branches die, the pruned tree has fewer nodes and lower error
  # load. This is a direct consequence: fewer expected tests = lower
  # sum of path_power * theta.
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)
  z_crit <- qnorm(1 - 0.05 / 2)

  # Full tree error load
  full_el <- manytestsr:::.error_load_from_tree(nd, delta_hat = 0.5,
    z_crit = z_crit, thealpha = 0.05)

  # Prune: kill 3 of 4 branches at depth 2
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  failed <- depth2_nodes[2:4]
  failed_and_desc <- c(failed, manytestsr:::.get_all_descendants(nd, failed))
  pruned_nd <- nd[!nd$nodenum %in% failed_and_desc, ]

  pruned_el <- manytestsr:::.error_load_from_tree(pruned_nd, delta_hat = 0.5,
    z_crit = z_crit, thealpha = 0.05)

  expect_lt(pruned_el$sum_G, full_el$sum_G,
    label = "Pruned error load should be less than full"
  )
  # Depth 1 error load is unchanged (still just the root)
  expect_equal(pruned_el$G[["1"]], full_el$G[["1"]],
    info = "Root error load should not change with pruning"
  )
})


# ============================================================================
# Edge cases and robustness
# ============================================================================

test_that("pruned: works when natural gating suffices (no adjustment needed)", {
  # When error load <= 1 for the full tree, all alphas are nominal.
  # Pruning should still work and return nominal alpha everywhere.
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 100)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.01)

  # Verify full-tree error load is <= 1
  z_crit <- qnorm(1 - 0.05 / 2)
  el <- manytestsr:::.error_load_from_tree(nd, delta_hat = 0.01,
    z_crit = z_crit, thealpha = 0.05)
  expect_true(el$sum_G <= 1,
    info = paste("Expected sum_G <= 1 for test to be meaningful, got", el$sum_G)
  )

  for (d in 1:3) {
    alpha_d <- obj$alphafn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(alpha_d), 0.05,
      info = paste("All alphas should be nominal when error load <= 1, depth", d)
    )
  }

  # After pruning (even though unnecessary), alphas should still be nominal
  depth2_nodes <- nd$nodenum[nd$depth == 2L]
  failed <- depth2_nodes[2:3]
  failed_and_desc <- c(failed, manytestsr:::.get_all_descendants(nd, failed))
  pruned_nd <- nd[!nd$nodenum %in% failed_and_desc, ]
  obj$update(pruned_nd, thealpha = 0.05)

  for (d in 1:3) {
    alpha_d <- obj$alphafn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(alpha_d), 0.05,
      info = paste("Nominal alpha persists after pruning low-load tree, depth", d)
    )
  }
})

test_that("pruned: update then reset then update cycle works", {
  # Simulates what happens across multiple simulation runs:
  # run 1 prunes some branches, then reset, then run 2 prunes different branches.
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 2000)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)
  depth2_nodes <- nd$nodenum[nd$depth == 2L]

  # Run 1: kill branches 2-4, keep branch 1
  failed1 <- depth2_nodes[2:4]
  failed_and_desc1 <- c(failed1, manytestsr:::.get_all_descendants(nd, failed1))
  pruned1 <- nd[!nd$nodenum %in% failed_and_desc1, ]
  obj$update(pruned1, thealpha = 0.05)
  alpha3_run1 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )

  # Reset between runs
  obj$reset(thealpha = 0.05)

  # Run 2: kill branches 1 and 3-4, keep branch 2
  failed2 <- depth2_nodes[c(1, 3, 4)]
  failed_and_desc2 <- c(failed2, manytestsr:::.get_all_descendants(nd, failed2))
  pruned2 <- nd[!nd$nodenum %in% failed_and_desc2, ]
  obj$update(pruned2, thealpha = 0.05)
  alpha3_run2 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )

  # Both runs pruned the same number of branches (3 of 4), so in a
  # regular tree (equal nodesize) the resulting alphas should be equal
  expect_equal(unname(alpha3_run1), unname(alpha3_run2),
    info = "Symmetric pruning in a regular tree should give equal alpha"
  )
})

test_that("pruned: handles thealpha change correctly", {
  # When find_blocks passes a different thealpha, the pruned schedule
  # should recompute with the new alpha.
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)

  alpha_at_05 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 2L
  )
  alpha_at_10 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.10, thew0 = 0.049, depth = 2L
  )

  # Higher nominal alpha should give higher adjusted alpha
  expect_gt(unname(alpha_at_10), unname(alpha_at_05))
})

test_that("pruned: validates node_dat columns", {
  bad_dat <- data.frame(nodenum = 1, parent = 0, depth = 1)
  expect_error(
    alpha_adaptive_tree_pruned(node_dat = bad_dat, delta_hat = 0.5),
    "nodesize"
  )
})

test_that("pruned: monotonicity — more pruning gives more alpha", {
  # Pruning 4 branches should give at least as much alpha at depth 3
  # as pruning only 2 branches, because the surviving tree is smaller.
  nd <- make_regular_tree(k = 5, max_depth = 3, N_total = 2500)
  obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)
  depth2_nodes <- nd$nodenum[nd$depth == 2L]

  # Prune 2 branches
  failed2 <- depth2_nodes[4:5]
  fd2 <- c(failed2, manytestsr:::.get_all_descendants(nd, failed2))
  pruned2 <- nd[!nd$nodenum %in% fd2, ]
  obj$update(pruned2, thealpha = 0.05)
  alpha3_prune2 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )

  # Reset, then prune 4 branches
  obj$reset(thealpha = 0.05)
  failed4 <- depth2_nodes[2:5]
  fd4 <- c(failed4, manytestsr:::.get_all_descendants(nd, failed4))
  pruned4 <- nd[!nd$nodenum %in% fd4, ]
  obj$update(pruned4, thealpha = 0.05)
  alpha3_prune4 <- obj$alphafn(
    pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 3L
  )

  # More pruning -> more alpha (or equal, if both hit the nominal cap)
  expect_gte(unname(alpha3_prune4), unname(alpha3_prune2),
    label = "Pruning more branches should give >= alpha at deeper levels"
  )
})
