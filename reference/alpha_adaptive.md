# Adaptive Alpha Adjustment Based on Power Decay

Factory function that creates an alpha adjustment function for use with
[`find_blocks`](find_blocks.md). The returned function adjusts
significance levels at each tree depth based on estimated power decay
(Algorithm 1 from Appendix B of the supplement).

## Usage

``` r
alpha_adaptive(k, delta_hat, N_total, max_depth = 20L)
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
once per unique value of `thealpha` and reused on subsequent calls. When
the error load is at most 1 (natural gating suffices), nominal alpha is
returned at every level.

## Examples

``` r
# Create an adaptive alpha function for a 4-ary tree
my_alpha <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)

# Use with find_blocks
# find_blocks(idat, bdat, ..., alphafn = my_alpha)

# Inspect the alpha schedule it will use
compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000)
#>            1            2            3            4            5            6 
#> 0.0500000000 0.0125000000 0.0031250000 0.0007997540 0.0003946938 0.0005959012 
#>            7            8            9           10           11           12 
#> 0.0020881421 0.0120383281 0.0500000000 0.0500000000 0.0500000000 0.0500000000 
#>           13           14           15           16           17           18 
#> 0.0500000000 0.0500000000 0.0500000000 0.0500000000 0.0500000000 0.0500000000 
#>           19           20 
#> 0.0500000000 0.0500000000 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>            1            2            3            4            5            6 
#> 1.000000e+00 4.000000e+00 1.562981e+01 3.167012e+01 2.097663e+01 5.986183e+00 
#>            7            8            9           10           11           12 
#> 1.038350e+00 1.376706e-01 1.587883e-02 1.706042e-03 1.768565e-04 1.800728e-05 
#>           13           14           15           16           17           18 
#> 1.817040e-06 1.825254e-07 1.829376e-08 1.831441e-09 1.832474e-10 1.832991e-11 
#>           19           20 
#> 1.833249e-12 1.833378e-13 
#> 
#> attr(,"error_load")$sum_G
#> [1] 80.45654
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>          1          2          3          4          5          6          7 
#> 1.00000000 1.00000000 0.97686287 0.50656612 0.16558692 0.07134347 0.04336445 
#>          8          9         10         11         12         13         14 
#> 0.03314649 0.02883482 0.02686032 0.02591620 0.02545465 0.02522646 0.02511302 
#>         15         16         17         18         19         20 
#> 0.02505646 0.02502821 0.02501410 0.02500705 0.02500353 0.02500176 
#> 
#> attr(,"error_load")$critical_level
#> [1] 5
#> 
#> attr(,"error_load")$n_by_level
#>            1            2            3            4            5            6 
#> 1.000000e+03 2.500000e+02 6.250000e+01 1.562500e+01 3.906250e+00 9.765625e-01 
#>            7            8            9           10           11           12 
#> 2.441406e-01 6.103516e-02 1.525879e-02 3.814697e-03 9.536743e-04 2.384186e-04 
#>           13           14           15           16           17           18 
#> 5.960464e-05 1.490116e-05 3.725290e-06 9.313226e-07 2.328306e-07 5.820766e-08 
#>           19           20 
#> 1.455192e-08 3.637979e-09 
#> 
```
