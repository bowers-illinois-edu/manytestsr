# Adaptive Alpha with Branch Pruning

Factory function that creates an adaptive alpha adjustment system
supporting branch pruning. Unlike
[`alpha_adaptive_tree`](alpha_adaptive_tree.md), which pre-computes a
fixed schedule, this version can recompute the schedule on a pruned
subtree after each depth — giving more alpha to surviving branches when
dead branches are removed.

## Usage

``` r
alpha_adaptive_tree_pruned(node_dat, delta_hat, max_depth = NULL)
```

## Arguments

- node_dat:

  A data.frame or data.table with columns `nodenum`, `parent`, `depth`,
  and `nodesize`. Typically extracted from a
  [`find_blocks`](find_blocks.md) result. The root node must have
  `parent = 0` and `depth = 1`.

- delta_hat:

  Estimated standardized effect size (e.g., Cohen's d). Conservative
  (larger) values produce more stringent adjustment, which preserves the
  FWER guarantee.

- max_depth:

  Maximum depth to compute. Defaults to the maximum depth present in
  `node_dat`.

## Value

A list with three components:

- `alphafn`:

  A function with the standard alphafn interface:
  `function(pval, batch, nodesize, thealpha, thew0, depth)`. Looks up
  the current alpha schedule by depth.

- `update`:

  A function `function(pruned_node_dat, thealpha)` that recomputes the
  alpha schedule on the given (pruned) tree. Called by `find_blocks`
  after each depth's testable decisions.

- `reset`:

  A function `function(thealpha)` that restores the alpha schedule to
  the full (unpruned) tree. Called by `find_blocks` at the start of each
  run to ensure independence across simulation iterations.

## Details

The FWER guarantee follows from the same telescoping-sum argument as
[`alpha_adaptive_tree`](alpha_adaptive_tree.md), applied to the
surviving subtree at each depth. Pruning decisions at depth \\d\\ depend
only on tests at depths \\1, \ldots, d\\, which are independent of tests
at deeper levels (under data splitting or independent permutation
tests), so the conditional FWER bound holds.

When `find_blocks` detects a list-valued `alphafn`, it extracts these
three components and calls `reset` at the start of each run and `update`
after each depth.

## Examples

``` r
nd <- data.frame(
  nodenum  = 1:7,
  parent   = c(0, 1, 1, 2, 2, 3, 3),
  depth    = c(1, 2, 2, 3, 3, 3, 3),
  nodesize = c(500, 250, 250, 125, 125, 100, 150)
)
obj <- alpha_adaptive_tree_pruned(node_dat = nd, delta_hat = 0.5)
# obj$alphafn — pass to find_blocks
# obj$update(pruned_nd, 0.05) — recompute on surviving tree
# obj$reset(0.05) — restore full-tree schedule
```
