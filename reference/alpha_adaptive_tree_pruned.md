# Adaptive Alpha with Branch Pruning

Factory function that creates an adaptive alpha adjustment system
supporting branch pruning. Unlike
[`alpha_adaptive_tree`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree.md),
which pre-computes a fixed schedule, this version can recompute the
schedule on a pruned subtree after each depth — giving more alpha to
surviving branches when dead branches are removed.

## Usage

``` r
alpha_adaptive_tree_pruned(
  node_dat,
  delta_hat,
  max_depth = NULL,
  budget_weights = NULL,
  budget_total = 1,
  spending_fraction = 0.5,
  switching = FALSE
)
```

## Arguments

- node_dat:

  A data.frame or data.table with columns `nodenum`, `parent`, `depth`,
  and `nodesize`. Typically extracted from a
  [`find_blocks`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)
  result. The root node must have `parent = 0` and `depth = 1`.

- delta_hat:

  Estimated standardized effect size (e.g., Cohen's d). Conservative
  (larger) values produce more stringent adjustment, which preserves the
  FWER guarantee.

- max_depth:

  Maximum depth to compute. Defaults to the maximum depth present in
  `node_dat`.

- budget_weights:

  Controls depth-wise budget allocation. Accepts the same values as
  [`compute_adaptive_alphas_tree`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas_tree.md),
  plus `"remaining"`: a sequential spending process where a fraction
  `spending_fraction` of the remaining budget is spent at each depth.
  See Details.

- budget_total:

  Initial error budget (default 1.0). The constraint \\\sum w\_\ell
  \le\\ `budget_total` guarantees FWER control.

- spending_fraction:

  Fraction of remaining budget to spend at each depth when
  `budget_weights = "remaining"` (default 0.5). At depth \\\ell\\, the
  weight is \\w\_\ell = f \times B\_\ell\\ where \\f\\ is the spending
  fraction and \\B\_\ell\\ is the remaining budget.

- switching:

  Logical (default `FALSE`). When `TRUE`, implements the switching
  corollary (Corollary B.1): after each update, if the remaining pruned
  error load fits within the remaining budget, all deeper depths revert
  to nominal alpha.

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

The FWER guarantee follows from Theorem B.5 in the supplement:
predictable budget weights with data-dependent denominators. The weights
are "predictable" because \\w\_\ell\\ depends only on the testing
history through depth \\\ell - 1\\, not on depth-\\\ell\\ outcomes. The
union bound across depths gives FWER \\\le \alpha\\ whenever \\\sum
w\_\ell \le 1\\.

The **switching corollary** (when `switching = TRUE`): after pruning
narrows the surviving tree, if the remaining error load \\\sum\_{\ell
\ge s} D\_\ell \le B_s\\ (remaining budget), then \\\alpha\_\ell =
\alpha\\ for all \\\ell \ge s\\. This works by setting \\w\_\ell =
D\_\ell\\, so the \\D\_\ell\\ in numerator and denominator cancel,
leaving nominal alpha.

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
# obj$alphafn -- pass to find_blocks
# obj$update(pruned_nd, 0.05) -- recompute on surviving tree
# obj$reset(0.05) -- restore full-tree schedule
```
