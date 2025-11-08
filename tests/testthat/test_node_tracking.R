# Tests for node tracking and ancestry functions

## Setup for tests
interactive <- FALSE
if (interactive) {
  library(testthat)
  local_edition(3)
  library(here)
  library(data.table)
  library(devtools)
  load_all()
}

# Access internal functions via namespace
create_node_tracker <- manytestsr:::create_node_tracker
add_nodes_to_tracker <- manytestsr:::add_nodes_to_tracker
get_parent_from_tracker <- manytestsr:::get_parent_from_tracker
build_numeric_ancestry <- manytestsr:::build_numeric_ancestry

# Helper function to add nodes and update tracker
add_nodes_helper <- function(tracker, parent_ids, depth, n_children) {
  result <- add_nodes_to_tracker(tracker, parent_ids, depth, n_children)
  result$tracker
}

# ============================================================================
# BASIC NODE TRACKER TESTS
# ============================================================================

test_that("create_node_tracker creates valid initial structure", {
  tracker <- create_node_tracker()

  # Should be a list with tracker and next_id components
  expect_type(tracker, "list")
  expect_true("tracker" %in% names(tracker))
  expect_true("next_id" %in% names(tracker))

  # Tracker should be a data.table
  expect_s3_class(tracker$tracker, "data.table")

  # Should have correct initial structure
  expect_equal(nrow(tracker$tracker), 1)
  expect_equal(names(tracker$tracker), c("node_id", "parent_id", "depth"))

  # Root node should have correct values
  expect_equal(tracker$tracker$node_id, 1L)
  expect_equal(tracker$tracker$parent_id, 0L)
  expect_equal(tracker$tracker$depth, 1L)

  # Next ID should be 2
  expect_equal(tracker$next_id, 2L)
})

test_that("add_nodes_to_tracker works correctly", {
  tracker <- create_node_tracker()

  # Add 2 children to root node (node 1) at depth 2
  result <- add_nodes_to_tracker(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)
  tracker <- result$tracker # Update tracker reference

  # Should now have 3 nodes total
  expect_equal(nrow(tracker$tracker), 3)
  expect_equal(tracker$next_id, 4L) # Next available ID should be 4

  # Check the new nodes
  new_nodes <- tracker$tracker[node_id %in% c(2L, 3L)]
  expect_equal(nrow(new_nodes), 2)
  expect_equal(new_nodes$parent_id, c(1L, 1L))
  expect_equal(new_nodes$depth, c(2L, 2L))
  expect_equal(new_nodes$node_id, c(2L, 3L))

  # Check return value includes new IDs
  expect_equal(result$new_ids, c(2L, 3L))
})

test_that("add_nodes_to_tracker works with multiple parents", {
  tracker <- create_node_tracker()

  # Add 2 children to root
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)

  # Add children to both of the new nodes
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 3L), depth = 3L, n_children = 4L)

  # Should now have 7 nodes total (1 root + 2 at depth 2 + 4 at depth 3)
  expect_equal(nrow(tracker$tracker), 7)
  expect_equal(tracker$next_id, 8L)

  # Check depth 3 nodes
  depth3_nodes <- tracker$tracker[depth == 3L]
  expect_equal(nrow(depth3_nodes), 4)
  expect_equal(depth3_nodes$node_id, c(4L, 5L, 6L, 7L))
  expect_equal(depth3_nodes$parent_id, c(2L, 3L, 2L, 3L)) # Round-robin assignment
})

# ============================================================================
# get_parent_from_tracker TESTS
# ============================================================================

test_that("get_parent_from_tracker returns correct parents", {
  tracker <- create_node_tracker()
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)
  tracker <- add_nodes_helper(tracker, parent_ids = 2L, depth = 3L, n_children = 2L)

  # Test single node lookups
  expect_equal(get_parent_from_tracker(tracker, 1L), 0L) # Root has no parent
  expect_equal(get_parent_from_tracker(tracker, 2L), 1L)
  expect_equal(get_parent_from_tracker(tracker, 3L), 1L)
  expect_equal(get_parent_from_tracker(tracker, 4L), 2L)
  expect_equal(get_parent_from_tracker(tracker, 5L), 2L)

  # Test multiple node lookups
  parents <- get_parent_from_tracker(tracker, c(2L, 3L, 4L, 5L))
  expect_equal(parents, c(1L, 1L, 2L, 2L))

  # Test non-existent node
  expect_equal(get_parent_from_tracker(tracker, 999L), 0L)
})

test_that("get_parent_from_tracker handles edge cases", {
  tracker <- create_node_tracker()

  # Empty vector - sapply returns list for empty input
  expect_equal(length(get_parent_from_tracker(tracker, integer(0))), 0)

  # Root node
  expect_equal(get_parent_from_tracker(tracker, 1L), 0L)

  # Non-existent nodes
  expect_equal(get_parent_from_tracker(tracker, c(100L, 200L)), c(0L, 0L))
})

# ============================================================================
# build_numeric_ancestry TESTS
# ============================================================================

test_that("build_numeric_ancestry works for simple tree", {
  tracker <- create_node_tracker()
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)
  tracker <- add_nodes_helper(tracker, parent_ids = 2L, depth = 3L, n_children = 2L)

  # Test root node - should return just itself
  expect_equal(build_numeric_ancestry(tracker, 1L), c(1L))

  # Test depth 2 nodes - should return path from root
  expect_equal(build_numeric_ancestry(tracker, 2L), c(1L, 2L))
  expect_equal(build_numeric_ancestry(tracker, 3L), c(1L, 3L))

  # Test depth 3 nodes - should return full path
  expect_equal(build_numeric_ancestry(tracker, 4L), c(1L, 2L, 4L))
  expect_equal(build_numeric_ancestry(tracker, 5L), c(1L, 2L, 5L))
})

test_that("build_numeric_ancestry works for deep tree", {
  tracker <- create_node_tracker()

  # Build a linear tree: 1 -> 2 -> 3 -> 4 -> 5
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 1L) # Node 2
  tracker <- add_nodes_helper(tracker, parent_ids = 2L, depth = 3L, n_children = 1L) # Node 3
  tracker <- add_nodes_helper(tracker, parent_ids = 3L, depth = 4L, n_children = 1L) # Node 4
  tracker <- add_nodes_helper(tracker, parent_ids = 4L, depth = 5L, n_children = 1L) # Node 5

  # Test deep ancestry
  expect_equal(build_numeric_ancestry(tracker, 5L), c(1L, 2L, 3L, 4L, 5L))
  expect_equal(build_numeric_ancestry(tracker, 4L), c(1L, 2L, 3L, 4L))
  expect_equal(build_numeric_ancestry(tracker, 3L), c(1L, 2L, 3L))
})

test_that("build_numeric_ancestry works for complex tree", {
  # Create a more complex tree structure
  tracker <- create_node_tracker()

  # Level 1: Root (1)
  # Level 2: Add 3 children (2, 3, 4)
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 3L)

  # Level 3: Add children to node 2 and node 4 (round-robin assignment)
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 4L), depth = 3L, n_children = 4L) # Nodes 5,6,7,8

  # Level 4: Add children to node 6
  tracker <- add_nodes_helper(tracker, parent_ids = 6L, depth = 4L, n_children = 2L) # Nodes 9,10

  # Test various paths - check what the actual structure is first
  expect_equal(build_numeric_ancestry(tracker, 1L), c(1L))
  expect_equal(build_numeric_ancestry(tracker, 3L), c(1L, 3L)) # Direct child of root

  # With round-robin assignment, node 6 should be child of node 4, not 2
  expect_equal(build_numeric_ancestry(tracker, 6L), c(1L, 4L, 6L))
  expect_equal(build_numeric_ancestry(tracker, 10L), c(1L, 4L, 6L, 10L)) # Deep path

  # Verify the tree structure is as expected
  expect_equal(nrow(tracker$tracker), 10) # Total nodes
  expect_equal(tracker$next_id, 11L) # Next available ID
})

test_that("build_numeric_ancestry handles edge cases", {
  tracker <- create_node_tracker()

  # Test root node
  expect_equal(build_numeric_ancestry(tracker, 1L), c(1L))

  # Test with character input (should be converted to integer)
  expect_equal(build_numeric_ancestry(tracker, "1"), c(1L))

  # Add some nodes for more tests
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)

  # Test numeric input as double
  expect_equal(build_numeric_ancestry(tracker, 2.0), c(1L, 2L))
})

# ============================================================================
# INTEGRATION TESTS WITH REALISTIC TREE STRUCTURES
# ============================================================================

test_that("node tracking works for binary tree structure", {
  # Create a balanced binary tree
  tracker <- create_node_tracker()

  # Level 2: 2 children
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)

  # Level 3: Each level 2 node gets 2 children
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 3L), depth = 3L, n_children = 4L)

  # Level 4: Each level 3 node gets 2 children
  tracker <- add_nodes_helper(tracker, parent_ids = c(4L, 5L, 6L, 7L), depth = 4L, n_children = 8L)

  # Test structure
  expect_equal(nrow(tracker$tracker), 15) # 1 + 2 + 4 + 8 = 15 nodes

  # Test ancestry for leaf nodes - with round-robin, paths will be different
  expect_equal(build_numeric_ancestry(tracker, 8L), c(1L, 2L, 4L, 8L)) # Left-most path
  expect_equal(build_numeric_ancestry(tracker, 15L), c(1L, 3L, 7L, 15L)) # Right-most path
  expect_equal(build_numeric_ancestry(tracker, 11L), c(1L, 3L, 7L, 11L)) # Middle path (corrected)

  # Test that all paths start with root
  for (node_id in 8L:15L) {
    path <- build_numeric_ancestry(tracker, node_id)
    expect_equal(path[1], 1L, info = paste("Node", node_id, "path should start with root"))
    expect_equal(length(path), 4L, info = paste("Node", node_id, "should have depth 4"))
  }
})

test_that("node tracking matches expected tree patterns", {
  # Test that our tracking system produces the expected hierarchical structure
  tracker <- create_node_tracker()

  # Build tree similar to what might be created by find_blocks
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 3L) # Nodes 2,3,4
  tracker <- add_nodes_helper(tracker, parent_ids = 2L, depth = 3L, n_children = 2L) # Nodes 5,6
  tracker <- add_nodes_helper(tracker, parent_ids = 4L, depth = 3L, n_children = 2L) # Nodes 7,8

  # Verify parent-child relationships
  children_of_1 <- tracker$tracker[parent_id == 1L, node_id]
  children_of_2 <- tracker$tracker[parent_id == 2L, node_id]
  children_of_4 <- tracker$tracker[parent_id == 4L, node_id]

  expect_equal(sort(children_of_1), c(2L, 3L, 4L))
  expect_equal(sort(children_of_2), c(5L, 6L))
  expect_equal(sort(children_of_4), c(7L, 8L))

  # Verify depths
  expect_equal(tracker$tracker[node_id %in% c(2L, 3L, 4L), unique(depth)], 2L)
  expect_equal(tracker$tracker[node_id %in% c(5L, 6L, 7L, 8L), unique(depth)], 3L)

  # Test ancestry paths reflect the structure
  expect_equal(build_numeric_ancestry(tracker, 5L), c(1L, 2L, 5L))
  expect_equal(build_numeric_ancestry(tracker, 8L), c(1L, 4L, 8L))

  # Node 3 is a leaf (has no children)
  expect_equal(build_numeric_ancestry(tracker, 3L), c(1L, 3L))
})

# ============================================================================
# PROPERTY-BASED TESTS
# ============================================================================

test_that("build_numeric_ancestry paths always start with root", {
  tracker <- create_node_tracker()
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 5L)
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 4L, 6L), depth = 3L, n_children = 6L)

  # Test that all ancestry paths start with 1 (root)
  all_nodes <- tracker$tracker$node_id
  for (node_id in all_nodes) {
    path <- build_numeric_ancestry(tracker, node_id)
    expect_equal(path[1], 1L, info = paste("Node", node_id, "ancestry should start with root"))
    expect_true(length(path) >= 1, info = paste("Node", node_id, "should have non-empty path"))
    expect_equal(path[length(path)], node_id, info = paste("Path should end with the node itself"))
  }
})

test_that("build_numeric_ancestry paths have consistent lengths with depth", {
  tracker <- create_node_tracker()
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 3L), depth = 3L, n_children = 4L)

  # Path length should equal the node's depth
  for (i in 1:nrow(tracker$tracker)) {
    node_row <- tracker$tracker[i, ]
    path <- build_numeric_ancestry(tracker, node_row$node_id)
    expect_equal(length(path), node_row$depth,
      info = paste("Node", node_row$node_id, "at depth", node_row$depth)
    )
  }
})

test_that("ancestry paths are consistent with parent relationships", {
  tracker <- create_node_tracker()
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 3L)
  tracker <- add_nodes_helper(tracker, parent_ids = c(2L, 3L), depth = 3L, n_children = 4L)

  # For every non-root node, its ancestry should be parent's ancestry + itself
  non_root_nodes <- tracker$tracker[node_id != 1L]
  for (i in 1:nrow(non_root_nodes)) {
    node_row <- non_root_nodes[i, ]
    node_path <- build_numeric_ancestry(tracker, node_row$node_id)
    parent_path <- build_numeric_ancestry(tracker, node_row$parent_id)

    # Node path should be parent path + node itself
    expected_path <- c(parent_path, node_row$node_id)
    expect_equal(node_path, expected_path,
      info = paste("Node", node_row$node_id, "ancestry inconsistent with parent")
    )
  }
})

# ============================================================================
# REALISTIC TREE TRACING TESTS
# ============================================================================

test_that("tree tracing works for find_blocks-like scenarios", {
  # Simulate a tree that might be created by find_blocks splitting
  tracker <- create_node_tracker()

  # Initial split: root -> 2 children (like left/right split)
  tracker <- add_nodes_helper(tracker, parent_ids = 1L, depth = 2L, n_children = 2L)

  # Second level splits: left child splits further
  tracker <- add_nodes_helper(tracker, parent_ids = 2L, depth = 3L, n_children = 2L) # Nodes 4,5

  # Third level: one more split on node 4
  tracker <- add_nodes_helper(tracker, parent_ids = 4L, depth = 4L, n_children = 2L) # Nodes 6,7

  # Test complete path tracing for deepest nodes
  path_6 <- build_numeric_ancestry(tracker, 6L)
  path_7 <- build_numeric_ancestry(tracker, 7L)

  expect_equal(path_6, c(1L, 2L, 4L, 6L))
  expect_equal(path_7, c(1L, 2L, 4L, 7L))

  # Test intermediate paths
  expect_equal(build_numeric_ancestry(tracker, 4L), c(1L, 2L, 4L))
  expect_equal(build_numeric_ancestry(tracker, 5L), c(1L, 2L, 5L))
  expect_equal(build_numeric_ancestry(tracker, 3L), c(1L, 3L)) # Right branch, no further splits

  # Verify tree structure integrity
  expect_equal(nrow(tracker$tracker), 7) # 1 + 2 + 2 + 2 = 7 nodes

  # Verify depth relationships
  depths <- tracker$tracker[, .(max_depth = max(depth)), by = node_id]
  root_depth <- depths[node_id == 1L, max_depth]
  deepest_depth <- max(depths$max_depth)

  expect_equal(root_depth, 1L)
  expect_equal(deepest_depth, 4L)
})
