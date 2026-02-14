# Tests for tree-based adaptive alpha (compute_adaptive_alphas_tree / alpha_adaptive_tree)
#
# These functions extend the parametric alpha_adaptive framework to handle
# irregular trees — where branching factor and sample sizes vary across nodes.
# The core statistical guarantee is the same: adjust alpha at each depth to
# compensate for the expected number of tests conducted at that depth.
#
# The denominator at depth d is sum(path_power) over all nodes at depth d,
# where path_power = product of ancestor thetas. This counts the expected
# number of tests reaching depth d. For a regular k-ary tree, this reduces
# to k^{d-1} * prod theta, matching the parametric formula.

# ============================================================================
# Helper: create node_dat for a regular k-ary tree
# ============================================================================

make_regular_tree <- function(k, max_depth, N_total) {
  # Builds a balanced k-ary tree where every node at depth d has

  # nodesize = N_total / k^{d-1} and exactly k children.
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
# Helper: DPP-like irregular tree
# ============================================================================

make_dpp_tree <- function() {
  # Mimics the DPP design: Root → 5 colleges → 12 cohorts → 16 leaves
  # with varying sample sizes at each level.
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

  # 5 colleges with varying sizes
  college_sizes <- c(450, 400, 450, 400, 500)
  college_ids <- integer(length(college_sizes))
  for (i in seq_along(college_sizes)) {
    nodes <- rbind(nodes, data.frame(
      nodenum = id, parent = root_id, depth = 2L, nodesize = college_sizes[i]
    ))
    college_ids[i] <- id; id <- id + 1L
  }

  # Cohorts within colleges (irregular: 2-3 per college)
  cohort_parents <- c(
    rep(college_ids[1], 3),  # college 1 → 3 cohorts
    rep(college_ids[2], 2),  # college 2 → 2 cohorts
    rep(college_ids[3], 3),  # college 3 → 3 cohorts
    rep(college_ids[4], 2),  # college 4 → 2 cohorts
    rep(college_ids[5], 2)   # college 5 → 2 cohorts
  )
  cohort_sizes <- c(150, 150, 150,    # college 1
                    200, 200,          # college 2
                    150, 150, 150,     # college 3
                    200, 200,          # college 4
                    200, 300)          # college 5
  cohort_ids <- integer(length(cohort_sizes))
  for (i in seq_along(cohort_sizes)) {
    nodes <- rbind(nodes, data.frame(
      nodenum = id, parent = cohort_parents[i], depth = 3L, nodesize = cohort_sizes[i]
    ))
    cohort_ids[i] <- id; id <- id + 1L
  }

  # Leaves within cohorts (1-2 per cohort)
  leaf_parents <- c(
    rep(cohort_ids[1], 2), cohort_ids[2],  # college 1 cohorts
    cohort_ids[3], cohort_ids[4],          # college 1 + college 2
    rep(cohort_ids[5], 2),                 # college 2
    cohort_ids[6], cohort_ids[7],          # college 3
    rep(cohort_ids[8], 2),                 # college 3
    cohort_ids[9], cohort_ids[10],         # college 4
    cohort_ids[11], cohort_ids[12]         # college 5
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
# Tests for compute_adaptive_alphas_tree()
# ============================================================================

test_that("tree: root alpha is always nominal regardless of tree shape", {
  # The root is a single test — no multiplicity adjustment needed,
  # regardless of what the tree looks like below it.
  for (k in c(2, 4)) {
    nd <- make_regular_tree(k = k, max_depth = 3, N_total = 1000)
    alphas <- compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, thealpha = 0.05
    )
    expect_equal(alphas[[1]], 0.05,
      info = paste("regular tree k =", k)
    )
  }
  # Also for the irregular DPP tree
  dpp <- make_dpp_tree()
  alphas_dpp <- compute_adaptive_alphas_tree(
    node_dat = dpp, delta_hat = 0.5, thealpha = 0.05
  )
  expect_equal(alphas_dpp[[1]], 0.05, info = "DPP tree")
})

test_that("tree: matches parametric for regular k-ary trees", {
  # A balanced k-ary tree with equal splitting should produce the same
  # alpha schedule as the parametric formula. This is the key equivalence
  # test: the tree-based computation generalizes the parametric one.
  for (k in c(2, 3, 4)) {
    nd <- make_regular_tree(k = k, max_depth = 4, N_total = 1000)
    tree_alphas <- compute_adaptive_alphas_tree(
      node_dat = nd, delta_hat = 0.5, thealpha = 0.05
    )
    param_alphas <- compute_adaptive_alphas(
      k = k, delta_hat = 0.5, N_total = 1000,
      max_depth = 4, thealpha = 0.05
    )
    # Compare at each depth (tree only goes to max_depth of node_dat)
    for (d in seq_along(tree_alphas)) {
      expect_equal(tree_alphas[[d]], param_alphas[[d]],
        tolerance = 1e-10,
        info = paste("k =", k, "depth =", d)
      )
    }
  }
})

test_that("tree: adjusted alpha never exceeds nominal alpha", {
  # The min(thealpha, ...) guard ensures we never give a level MORE
  # alpha than nominal. This preserves FWER control.
  dpp <- make_dpp_tree()
  alphas <- compute_adaptive_alphas_tree(
    node_dat = dpp, delta_hat = 0.8, thealpha = 0.05
  )
  expect_true(all(unname(alphas) <= 0.05))
})

test_that("tree: all alphas are positive", {
  dpp <- make_dpp_tree()
  alphas <- compute_adaptive_alphas_tree(
    node_dat = dpp, delta_hat = 0.5, thealpha = 0.05
  )
  expect_true(all(alphas > 0))
})

test_that("tree: natural gating bypass when sum_G <= 1", {
  # When the total error load is at most 1, natural gating suffices
  # and all depths should get nominal alpha (no adjustment needed).
  # Use tiny delta_hat to ensure low power → low error load.
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 100)
  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.01, thealpha = 0.05
  )
  # Verify error load is actually <= 1
  el <- attr(alphas, "error_load")
  expect_true(el$sum_G <= 1,
    info = paste("Expected sum_G <= 1 but got", el$sum_G)
  )
  # All alphas should be nominal
  expect_true(all(unname(alphas) == 0.05))
})

test_that("tree: irregular tree differs from parametric approximation", {
  # The DPP tree is irregular — using k=3 as an approximation should
  # give different alphas than using the actual tree structure.
  dpp <- make_dpp_tree()
  tree_alphas <- compute_adaptive_alphas_tree(
    node_dat = dpp, delta_hat = 0.8, thealpha = 0.05
  )
  param_alphas <- compute_adaptive_alphas(
    k = 3, delta_hat = 0.8, N_total = 2200,
    max_depth = length(tree_alphas), thealpha = 0.05
  )
  # They should differ at depth 2 or deeper because the tree is irregular
  # (5 children at root, not 3; varying sizes throughout)
  diffs <- unname(tree_alphas) - unname(param_alphas[seq_along(tree_alphas)])
  expect_false(all(abs(diffs) < 1e-10),
    info = "Irregular tree should produce different alphas than k=3 parametric"
  )
})

test_that("tree: validates node_dat columns", {
  # Missing required columns should produce a clear error message.
  bad_dat <- data.frame(nodenum = 1, parent = 0, depth = 1)
  expect_error(
    compute_adaptive_alphas_tree(node_dat = bad_dat, delta_hat = 0.5),
    "nodesize"
  )
  bad_dat2 <- data.frame(nodenum = 1, nodesize = 100, depth = 1)
  expect_error(
    compute_adaptive_alphas_tree(node_dat = bad_dat2, delta_hat = 0.5),
    "parent"
  )
})

test_that("tree: has error_load attribute", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )
  el <- attr(alphas, "error_load")
  expect_true(!is.null(el))
  expect_true("sum_G" %in% names(el))
  expect_true("G" %in% names(el))
})

test_that("tree: deeper levels get stricter alpha when power is high", {
  # With large N and high effect, the procedure should tighten alpha
  # at deeper levels because many tests will be conducted.
  nd <- make_regular_tree(k = 4, max_depth = 3, N_total = 5000)
  alphas <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.8, thealpha = 0.05
  )
  # Each level should have equal or smaller alpha than the one above
  for (d in 2:length(alphas)) {
    expect_true(alphas[[d]] <= alphas[[d - 1]],
      info = paste("depth", d, "should be <= depth", d - 1)
    )
  }
})


# ============================================================================
# Tests for alpha_adaptive_tree() factory
# ============================================================================

test_that("tree factory: returns function with correct interface", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)
  expect_true(is.function(fn))
  # Must have the same formals as alpha_adaptive()
  args <- names(formals(fn))
  expect_true("pval" %in% args)
  expect_true("batch" %in% args)
  expect_true("nodesize" %in% args)
  expect_true("thealpha" %in% args)
  expect_true("thew0" %in% args)
  expect_true("depth" %in% args)
})

test_that("tree factory: consistent with compute_adaptive_alphas_tree", {
  nd <- make_regular_tree(k = 3, max_depth = 4, N_total = 1000)
  schedule <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)

  for (d in seq_along(schedule)) {
    result <- fn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(result), schedule[[d]], info = paste("depth", d))
  }
})

test_that("tree factory: handles mixed depths in single call", {
  nd <- make_regular_tree(k = 3, max_depth = 4, N_total = 1000)
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)
  schedule <- compute_adaptive_alphas_tree(
    node_dat = nd, delta_hat = 0.5, thealpha = 0.05
  )

  # Simulate an ancestry path with nodes at depths 1, 2, 3
  result <- fn(
    pval = c(0.01, 0.02, 0.03),
    batch = c(1, 2, 3),
    nodesize = c(1000, 333, 111),
    thealpha = 0.05,
    thew0 = 0.049,
    depth = c(1L, 2L, 3L)
  )
  expect_length(result, 3)
  expect_equal(unname(result[1]), schedule[[1]])
  expect_equal(unname(result[2]), schedule[[2]])
  expect_equal(unname(result[3]), schedule[[3]])
})

test_that("tree factory: caches correctly and recomputes on thealpha change", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)

  # Same thealpha → same result
  res1 <- fn(pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  res2 <- fn(pval = 0.99, batch = 1, nodesize = 100,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  expect_equal(res1, res2)

  # Different thealpha → different result
  res3 <- fn(pval = 0.01, batch = 1, nodesize = 100,
    thealpha = 0.10, thew0 = 0.049, depth = 2L)
  expect_false(isTRUE(all.equal(res1, res3)))
})

test_that("tree factory: ignores pval (alpha depends only on tree structure)", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)

  res_small_p <- fn(pval = 0.001, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  res_large_p <- fn(pval = 0.99, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  expect_equal(res_small_p, res_large_p)
})

test_that("tree factory: warns when depth is NULL", {
  nd <- make_regular_tree(k = 3, max_depth = 3, N_total = 1000)
  fn <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)
  expect_warning(
    fn(pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = NULL),
    "depth"
  )
})
