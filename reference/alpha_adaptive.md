# Adaptive Alpha Adjustment Based on Power Decay

Factory function that creates an alpha adjustment function for use with
[`find_blocks`](find_blocks.md). The returned function adjusts
significance levels at each tree depth based on estimated power decay
(Algorithm 1 from the adaptive alpha appendix).

## Usage

``` r
alpha_adaptive(k, delta_hat, N_total, tau = 0.1, max_depth = 20L)
```

## Arguments

- k:

  Branching factor. Either a scalar (constant k at all levels) or an
  integer vector of length `max_depth` where `k[ell]` is the branching
  factor at level `ell`.

- delta_hat:

  Estimated standardized effect size (e.g., Cohen's d). Conservative
  (larger) values produce more stringent adjustment, which preserves the
  FWER guarantee. Use an upper bound on the true effect size.

- N_total:

  Total sample size at the root level.

- tau:

  Cumulative power threshold (default 0.1). When cumulative power drops
  below `tau`, natural gating is deemed sufficient and nominal alpha is
  used.

- max_depth:

  Maximum depth to compute (default 20).

## Value

A function with signature
`function(pval, batch, nodesize, thealpha, thew0, depth)` conforming to
the `alphafn` interface used by [`find_blocks`](find_blocks.md).

## Details

The returned function uses the `depth` parameter (passed by
`find_blocks`) to look up the pre-computed alpha for each node's tree
depth. The `pval`, `batch`, `nodesize`, and `thew0` parameters are
accepted for interface compatibility but are not used — unlike online
FDR methods, the adaptive alpha depends only on tree structure, not on
observed p-values.

Results are cached internally: the vector of adjusted alphas is computed
once per unique value of `thealpha` and reused on subsequent calls.

## Examples

``` r
# Create an adaptive alpha function for a 4-ary tree
my_alpha <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)

# Use with find_blocks
# find_blocks(idat, bdat, ..., alphafn = my_alpha)

# Inspect the alpha schedule it will use
compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000)
#>           1           2           3           4           5           6 
#> 0.050000000 0.012500000 0.003125000 0.000799754 0.050000000 0.050000000 
#>           7           8           9          10          11          12 
#> 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 
#>          13          14          15          16          17          18 
#> 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 
#>          19          20 
#> 0.050000000 0.050000000 
```
