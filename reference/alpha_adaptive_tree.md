# Adaptive Alpha from an Actual Tree Structure

Factory function counterpart to [`alpha_adaptive`](alpha_adaptive.md)
for irregular trees. Creates an alpha adjustment function for use with
[`find_blocks`](find_blocks.md), using actual per-node sample sizes to
compute the alpha schedule instead of assuming a regular k-ary tree.

## Usage

``` r
alpha_adaptive_tree(node_dat, delta_hat, max_depth = NULL)
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

A function with signature
`function(pval, batch, nodesize, thealpha, thew0, depth)` conforming to
the `alphafn` interface used by [`find_blocks`](find_blocks.md).

## Details

The returned function behaves identically to the closure from
[`alpha_adaptive`](alpha_adaptive.md): it uses the `depth` parameter to
look up pre-computed alphas, ignoring `pval`, `batch`, and other
arguments. The difference is that the alpha schedule comes from
[`compute_adaptive_alphas_tree`](compute_adaptive_alphas_tree.md) rather
than the parametric formula, so it handles irregular branching
correctly.

Results are cached internally: the alpha schedule is computed once per
unique value of `thealpha` and reused on subsequent calls.

## Examples

``` r
# Build a tree template from the DPP design
nd <- data.frame(
  nodenum  = 1:7,
  parent   = c(0, 1, 1, 2, 2, 3, 3),
  depth    = c(1, 2, 2, 3, 3, 3, 3),
  nodesize = c(500, 250, 250, 125, 125, 100, 150)
)
my_alpha <- alpha_adaptive_tree(node_dat = nd, delta_hat = 0.5)

# Use with find_blocks
# find_blocks(idat, bdat, ..., alphafn = my_alpha)
```
