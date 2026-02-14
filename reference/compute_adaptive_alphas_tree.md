# Compute Adaptive Alpha Levels from an Actual Tree

Tree-mode counterpart to
[`compute_adaptive_alphas`](compute_adaptive_alphas.md). Instead of
assuming a regular k-ary tree with equal splits, this function takes the
actual tree structure (with per-node sample sizes) and computes
per-depth adjusted significance levels. This handles irregular trees
where branching factor and sample sizes vary across nodes.

## Usage

``` r
compute_adaptive_alphas_tree(
  node_dat,
  delta_hat,
  max_depth = NULL,
  thealpha = 0.05
)
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

- thealpha:

  Nominal significance level (default 0.05).

## Value

Named numeric vector of adjusted alpha levels, one per depth (1 through
`max_depth`). Names are depth levels as characters. Has attribute
`"error_load"` containing the
[`compute_error_load`](compute_error_load.md) result from the tree.

## Details

The algorithm mirrors
[`compute_adaptive_alphas`](compute_adaptive_alphas.md) but uses actual
per-node power instead of the parametric assumption:

1.  Compute per-node power \\\theta\\ and path power (product of
    ancestor thetas) via `compute_error_load`.

2.  If total error load \\\sum G\_\ell \le 1\\: natural gating suffices,
    return nominal `thealpha` at every depth.

3.  Otherwise, at depth \\d\\: \$\$\alpha_d = \min\left\\\alpha,\\
    \frac{\alpha}{\sum\_{\text{nodes at depth } d}
    \text{path\\power}(\text{node})} \right\\\$\$

The denominator at depth \\d\\ is the sum of path powers over all nodes
at that depth — i.e., the expected number of tests the procedure
conducts at depth \\d\\. For a regular k-ary tree with equal sample
sizes, this reduces to \\k^{d-1} \prod\_{j=1}^{d-1} \theta_j\\, matching
the parametric formula in
[`compute_adaptive_alphas`](compute_adaptive_alphas.md).

## Examples

``` r
# Build a small irregular tree
nd <- data.frame(
  nodenum  = 1:7,
  parent   = c(0, 1, 1, 2, 2, 3, 3),
  depth    = c(1, 2, 2, 3, 3, 3, 3),
  nodesize = c(500, 250, 250, 125, 125, 100, 150)
)
compute_adaptive_alphas_tree(node_dat = nd, delta_hat = 0.5)
#>      1      2      3 
#> 0.0500 0.0250 0.0125 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>        1        2        3 
#> 1.000000 2.000000 3.998518 
#> 
#> attr(,"error_load")$sum_G
#> [1] 6.998518
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>         1         2         3 
#> 1.0000000 1.0000000 0.9996296 
#> 
#> attr(,"error_load")$critical_level
#> [1] NA
#> 
#> attr(,"error_load")$n_by_level
#>   1   2   3 
#> 500 250 125 
#> 
#> attr(,"error_load")$node_detail
#>   nodenum parent depth nodesize     theta path_power
#> 1       1      0     1      500 1.0000000          1
#> 2       2      1     2      250 1.0000000          1
#> 3       3      1     2      250 1.0000000          1
#> 4       4      2     3      125 0.9998584          1
#> 5       5      2     3      125 0.9998584          1
#> 6       6      3     3      100 0.9988173          1
#> 7       7      3     3      150 0.9999843          1
#> 
```
