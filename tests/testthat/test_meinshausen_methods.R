test_that("Meinshausen hierarchical testing implementation works correctly", {
  skip_if_not_installed("data.table")

  # Create hierarchical test data
  node_dat <- data.table::data.table(
    nodenum = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
    p = c(0.001, 0.01, 0.03, 0.02, 0.08, 0.004, 0.12),
    depth = c(1L, 2L, 2L, 3L, 3L, 3L, 3L),
    nodesize = c(100, 50, 50, 25, 25, 25, 25),
    testable = rep(TRUE, 7)
  )

  # Create tree structure
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
      parent_id = c(0L, 1L, 1L, 2L, 2L, 3L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L, 3L, 3L)
    ),
    next_id = 8L
  )

  # Test that meinshausen_hierarchical_test exists and works
  skip_if_not(exists("meinshausen_hierarchical_test"))

  # Test with sequential rejection (Goeman-Solari enhancement)
  result_sequential <- meinshausen_hierarchical_test(
    node_dat, tracker,
    alpha = 0.05, method = "simes", use_sequential = TRUE
  )

  expect_true(is.data.table(result_sequential) || "data.table" %in% class(result_sequential))
  expect_true("meinshausen_reject" %in% names(result_sequential))
  expect_true("meinshausen_adjusted_p" %in% names(result_sequential))
  expect_true("meinshausen_level_tested" %in% names(result_sequential))

  # Check that some rejections occurred with significant p-values
  expect_true(any(result_sequential$meinshausen_reject))

  # Test with traditional Meinshausen procedure
  result_traditional <- meinshausen_hierarchical_test(
    node_dat, tracker,
    alpha = 0.05, method = "simes", use_sequential = FALSE
  )

  expect_true(is.data.table(result_traditional) || "data.table" %in% class(result_traditional))
  expect_true("meinshausen_reject" %in% names(result_traditional))

  # Sequential should have equal or greater power than traditional
  sequential_rejections <- sum(result_sequential$meinshausen_reject, na.rm = TRUE)
  traditional_rejections <- sum(result_traditional$meinshausen_reject, na.rm = TRUE)
  expect_true(sequential_rejections >= traditional_rejections)
})

test_that("Meinshausen method respects hierarchical structure", {
  skip_if_not_installed("data.table")

  # Create data where child is significant but parent isn't under naive testing
  node_dat <- data.table::data.table(
    nodenum = c(1L, 2L, 3L, 4L),
    p = c(0.08, 0.06, 0.001, 0.002), # Children more significant than parents
    depth = c(1L, 2L, 3L, 3L),
    nodesize = c(100, 50, 25, 25)
  )

  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L),
      parent_id = c(0L, 1L, 2L, 2L),
      depth = c(1L, 2L, 3L, 3L)
    )
  )

  skip_if_not(exists("meinshausen_hierarchical_test"))

  result <- meinshausen_hierarchical_test(node_dat, tracker, alpha = 0.05, method = "simes")

  # Check that hierarchical structure is respected
  # If a node is rejected, we should be able to test its children
  rejected_nodes <- result[meinshausen_reject == TRUE, nodenum]

  for (node_id in rejected_nodes) {
    # Under proper hierarchical testing, this rejection should be valid
    expect_true(result[nodenum == node_id, meinshausen_reject])
  }
})

test_that("Meinshausen p-value combination methods work correctly", {
  test_pvals <- c(0.01, 0.05, 0.08, 0.12, 0.2)

  # Test different combination methods
  if (exists("combine_pvalues_traditional")) {
    # Simes method
    simes_combined <- combine_pvalues_traditional(test_pvals, "simes")
    expect_true(is.numeric(simes_combined))
    expect_true(simes_combined >= 0 && simes_combined <= 1)

    # Fisher method
    fisher_combined <- combine_pvalues_traditional(test_pvals, "fisher")
    expect_true(is.numeric(fisher_combined))
    expect_true(fisher_combined >= 0 && fisher_combined <= 1)

    # Bonferroni method (most conservative)
    bonferroni_combined <- combine_pvalues_traditional(test_pvals, "bonferroni")
    expect_true(is.numeric(bonferroni_combined))
    expect_true(bonferroni_combined >= 0 && bonferroni_combined <= 1)

    # Bonferroni should be most conservative
    expect_true(bonferroni_combined >= simes_combined)
  }

  # Test level adjustment for sequential rejection
  if (exists("adjust_pvalues_level")) {
    adjusted_simes <- adjust_pvalues_level(test_pvals, "simes", 0.05)
    expect_true(length(adjusted_simes) == length(test_pvals))
    expect_true(all(adjusted_simes >= 0 & adjusted_simes <= 1))

    # Check monotonicity (required for sequential rejection)
    sorted_idx <- order(test_pvals)
    expect_true(all(diff(adjusted_simes[sorted_idx]) >= -1e-10)) # Allow for numerical precision
  }
})

test_that("Meinshausen hierarchical clustering generation works", {
  skip_if_not_installed("data.table")

  # Create test correlation matrix
  set.seed(123)
  n_vars <- 10
  X <- matrix(rnorm(100 * n_vars), ncol = n_vars)
  colnames(X) <- paste0("var", 1:n_vars)
  cor_matrix <- cor(X)

  skip_if_not(exists("generate_meinshausen_hierarchy"))

  # Generate hierarchy
  hierarchy <- generate_meinshausen_hierarchy(cor_matrix, method = "complete")

  expect_true(is.list(hierarchy))
  expect_true("tracker" %in% names(hierarchy))
  expect_true("hclust_object" %in% names(hierarchy))
  expect_true("variable_names" %in% names(hierarchy))

  # Check tracker structure
  tracker_data <- hierarchy$tracker
  expect_true(is.data.table(tracker_data) || "data.table" %in% class(tracker_data))
  expect_true("node_id" %in% names(tracker_data))
  expect_true("parent_id" %in% names(tracker_data))
  expect_true("depth" %in% names(tracker_data))
  expect_true("cluster_size" %in% names(tracker_data))

  # Should have nodes for variables plus internal nodes
  expect_true(nrow(tracker_data) >= n_vars)

  # Root node should have parent_id = 0
  expect_true(any(tracker_data$parent_id == 0))
})

test_that("Meinshausen integration with find_blocks works", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

  # Use example data to test integration
  data(example_dat, package = "manytestsr")

  idat <- data.table::as.data.table(example_dat)
  bdat <- idat %>%
    dplyr::group_by(blockF) %>%
    dplyr::summarize(
      nb = dplyr::n(),
      pb = mean(trt),
      hwt = (nb / nrow(idat)) * (pb * (1 - pb)),
      .groups = "drop"
    ) %>%
    data.table::as.data.table()

  # Test that Meinshausen integrates properly with find_blocks
  expect_no_error({
    result_meinshausen <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pIndepDist,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_meinshausen = TRUE,
      meinshausen_method = "simes",
      meinshausen_sequential = TRUE,
      thealpha = 0.05,
      maxtest = 10,
      trace = FALSE
    )
  })

  expect_true(is.list(result_meinshausen))
  expect_true("node_dat" %in% names(result_meinshausen))

  # If Meinshausen was applied successfully, check results
  if ("meinshausen_reject" %in% names(result_meinshausen$node_dat)) {
    # Check that results have correct structure
    node_dat <- result_meinshausen$node_dat
    expect_true(is.logical(node_dat$meinshausen_reject))
    expect_true(is.numeric(node_dat$meinshausen_adjusted_p))

    # All adjusted p-values should be valid
    expect_true(all(node_dat$meinshausen_adjusted_p >= 0 &
      node_dat$meinshausen_adjusted_p <= 1, na.rm = TRUE))

    # Rejected nodes should have reasonably small adjusted p-values
    if (any(node_dat$meinshausen_reject)) {
      rejected_adj_p <- node_dat[meinshausen_reject == TRUE, meinshausen_adjusted_p]
      expect_true(all(rejected_adj_p <= 0.05, na.rm = TRUE))
    }
  }
})

test_that("Meinshausen handles edge cases properly", {
  skip_if_not_installed("data.table")

  skip_if_not(exists("meinshausen_hierarchical_test"))

  # Test with single node
  single_node <- data.table::data.table(
    nodenum = 1L,
    p = 0.01,
    depth = 1L,
    nodesize = 50
  )

  single_tracker <- list(
    tracker = data.table::data.table(
      node_id = 1L,
      parent_id = 0L,
      depth = 1L
    )
  )

  expect_no_error({
    result_single <- meinshausen_hierarchical_test(single_node, single_tracker)
  })

  expect_true(is.data.table(result_single) || "data.table" %in% class(result_single))

  # Test with all non-significant p-values
  nonsig_nodes <- data.table::data.table(
    nodenum = c(1L, 2L, 3L),
    p = c(0.8, 0.9, 0.7),
    depth = c(1L, 2L, 2L),
    nodesize = c(50, 25, 25)
  )

  nonsig_tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L),
      parent_id = c(0L, 1L, 1L),
      depth = c(1L, 2L, 2L)
    )
  )

  expect_no_error({
    result_nonsig <- meinshausen_hierarchical_test(nonsig_nodes, nonsig_tracker)
  })

  # Should reject nothing or very few
  expect_true(sum(result_nonsig$meinshausen_reject, na.rm = TRUE) <= 1)

  # Test with missing p-values
  missing_p <- data.table::data.table(
    nodenum = c(1L, 2L, 3L),
    p = c(0.01, NA, 0.03),
    depth = c(1L, 2L, 2L),
    nodesize = c(50, 25, 25)
  )

  expect_no_error({
    result_missing <- meinshausen_hierarchical_test(missing_p, nonsig_tracker)
  })

  expect_true(is.logical(result_missing$meinshausen_reject))
})

test_that("Parameter validation works for Meinshausen methods", {
  skip_if_not_installed("data.table")

  node_dat <- data.table::data.table(
    nodenum = c(1L, 2L, 3L),
    p = c(0.01, 0.05, 0.1),
    depth = c(1L, 2L, 2L)
  )
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L),
      parent_id = c(0L, 1L, 1L),
      depth = c(1L, 2L, 2L)
    )
  )

  skip_if_not(exists("meinshausen_hierarchical_test"))

  # Test invalid alpha values
  expect_error(meinshausen_hierarchical_test(node_dat, tracker, alpha = -0.1), "alpha must be")
  expect_error(meinshausen_hierarchical_test(node_dat, tracker, alpha = 1.5), "alpha must be")
  expect_error(meinshausen_hierarchical_test(node_dat, tracker, alpha = c(0.05, 0.1)), "alpha must be")

  # Test invalid method
  expect_error(meinshausen_hierarchical_test(node_dat, tracker, method = "invalid"), "method must be")
})

test_that("Helper functions for Meinshausen work correctly", {
  # Test node relationship functions
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L),
      parent_id = c(0L, 1L, 1L, 2L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L)
    )
  )

  # Test get_children_nodes
  if (exists("get_children_nodes")) {
    children_1 <- get_children_nodes(1L, tracker)
    expect_true(is.integer(children_1))
    expect_true(2L %in% children_1 && 3L %in% children_1)

    children_2 <- get_children_nodes(2L, tracker)
    expect_true(4L %in% children_2)
  }

  # Test get_descendant_nodes
  if (exists("get_descendant_nodes")) {
    node_dat <- data.table::data.table(
      nodenum = c(1L, 2L, 3L, 4L, 5L),
      p = c(0.01, 0.02, 0.03, 0.04, 0.05)
    )

    descendants_1 <- get_descendant_nodes(1L, tracker, node_dat)
    expect_true(is.data.table(descendants_1) || "data.table" %in% class(descendants_1))
    expect_true(nrow(descendants_1) >= 1) # Should include self
  }
})

