## test_node_truth_propagation.R
##
## make_results_tree() labels block-level truth (`nonnull`) on LEAF nodes only,
## because blocks live at the leaves. Internal nodes were left NA, so the
## node-level FWER tally (node_any_false_rejection and the other node_* metrics)
## silently DROPPED any false rejection of an internal null hypothesis -- e.g.
## rejecting a whole null college or region. In structured top-down testing the
## procedure tests those internal hypotheses, so rejecting a true internal null IS
## a Type I error and must count toward the node-level FWER. (This is the bug that
## made the simulated Detroit Promise node-level FWER read as equal to the leaf
## FWER, ~0.04, when the true value with the internal Schoolcraft-college false
## rejection counted is ~0.08.)
##
## These tests pin the substantive requirement: truth is propagated up the tree
## (an internal node is non-null iff a descendant leaf is, null iff all its
## descendant leaves are known nulls), and a false rejection at an internal null
## node IS counted at the node level even when no leaf is falsely rejected.

library(data.table)

## A controlled 3-level tree:
##   root (1) -> college A (2, non-null) -> leaves 4 (effect), 5 (null)
##            -> college B (3, NULL)     -> leaves 6 (null),   7 (null)
## College B is a true null (all its blocks are null) yet is rejected at p = 0.01
## <= 0.05 -- a false rejection of an internal hypothesis -- while NO leaf is
## falsely rejected (leaves 5, 6, 7 all have p = 0.5 > 0.05). Leaf-only truth
## would miss B; propagated truth must catch it.
make_toy_truth_tree <- function() {
  node_dat <- data.table(
    nodenum = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
    parent  = c(0L, 1L, 1L, 2L, 2L, 3L, 3L),
    depth   = c(1L, 2L, 2L, 3L, 3L, 3L, 3L),
    p       = c(0.001, 0.001, 0.010, 0.020, 0.500, 0.500, 0.500),
    a       = rep(0.05, 7L)
  )
  ## Blocks live at the four leaves; only b1 (under college A) carries an effect.
  bdat <- data.table(
    blockF          = factor(c("b1", "b2", "b3", "b4")),
    nodenum_current = c(4L, 5L, 6L, 7L),
    trueblocks      = c(1, 0, 0, 0)
  )
  list(node_dat = node_dat, bdat = bdat)
}

run_toy <- function(truevar_name = "trueblocks") {
  toy <- make_toy_truth_tree()
  make_results_tree(
    orig_res = toy$bdat, block_id = "blockF",
    truevar_name = truevar_name, node_dat = toy$node_dat,
    return_what = c("nodes", "test_summary")
  )
}

test_that("internal-node truth is propagated up the tree", {
  out <- run_toy()
  nn  <- setNames(out$nodes$nonnull, as.character(out$nodes$name))
  expect_true(isTRUE(nn[["2"]]))          # college A contains an effect -> non-null
  expect_false(is.na(nn[["3"]]))          # college B is labeled, not left NA
  expect_false(isTRUE(nn[["3"]]))         # ... and it is a known null
  expect_true(isTRUE(nn[["1"]]))          # root contains an effect -> non-null
})

test_that("a false rejection of an internal null node counts toward node FWER", {
  out <- run_toy()
  ts  <- out$test_summary
  ## College B (internal null) is rejected at p = 0.01 <= 0.05: a node-level false
  ## rejection, with no leaf falsely rejected. The node tally must see it; the leaf
  ## tally must not. This is the discrepancy the bug hid.
  expect_true(isTRUE(ts$node_any_false_rejection))
  expect_false(isTRUE(ts$leaf_any_false_rejection))
  expect_equal(as.numeric(ts$node_num_false_rejections), 1)  # exactly college B
})

test_that("node_true_discoveries counts non-null nodes at every level", {
  out <- run_toy()
  ts  <- out$test_summary
  ## Non-null nodes rejected (p <= a): root (1), college A (2), leaf 4 -> 3.
  ## Leaf-level true discoveries: only leaf 4 among the leaves.
  expect_equal(as.numeric(ts$node_true_discoveries), 3)
  expect_equal(as.numeric(ts$leaf_true_discoveries), 1)
})

test_that("propagation invents no truth when there is none", {
  ## With no truevar, every node's nonnull stays NA and the node metrics are not
  ## computed; propagation must not manufacture FALSE/TRUE labels out of NA.
  out <- run_toy(truevar_name = NULL)
  expect_true(all(is.na(out$nodes$nonnull)))
})
