# Compute Adaptive Alpha Levels by Tree Depth

Implements Algorithm 1 from Appendix B of the supplement. First checks
whether natural gating suffices (total error load \\\le 1\\); if so,
returns nominal alpha at every level. Otherwise, computes adjusted
significance levels that compensate for the error load at each depth.

## Usage

``` r
compute_adaptive_alphas(
  k,
  delta_hat,
  N_total,
  max_depth = 20L,
  thealpha = 0.05
)
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

- thealpha:

  Nominal significance level (default 0.05).

## Value

Named numeric vector of adjusted alpha levels, one per depth (1 through
`max_depth`). Names are depth levels as characters. Has attribute
`"error_load"` containing the
[`compute_error_load`](compute_error_load.md) result, so the caller can
inspect whether adjustment was needed.

## Details

The function first calls [`compute_error_load`](compute_error_load.md)
to assess whether natural gating suffices. When \\\sum G\_\ell \le 1\\,
no adjustment is needed and nominal `thealpha` is returned at every
level.

When adjustment is needed, the formula at level \\\ell\\ is:
\$\$\alpha\_\ell^{adj} = \min\left\\\alpha,\\ \frac{\alpha}{k^{(\ell-1)}
\cdot \prod\_{j=1}^{\ell-1} \hat\theta_j} \right\\\$\$

The FWER guarantee (Theorem in the supplement) requires that power is
not underestimated (i.e., \\\hat\theta_j \geq \theta_j\\). In practice
this means using a conservatively large `delta_hat`.

## Examples

``` r
# Natural gating sufficient: all alphas = 0.05
compute_adaptive_alphas(k = 3, delta_hat = 0.2, N_total = 100,
                        max_depth = 4)
#>    1    2    3    4 
#> 0.05 0.05 0.05 0.05 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>          1          2          3          4 
#> 0.51596779 0.32557645 0.09567467 0.01653857 
#> 
#> attr(,"error_load")$sum_G
#> [1] 0.9537575
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] FALSE
#> 
#> attr(,"error_load")$thetas
#>          1          2          3          4 
#> 0.51596779 0.21033384 0.09795412 0.05762086 
#> 
#> attr(,"error_load")$critical_level
#> [1] 2
#> 
#> attr(,"error_load")$n_by_level
#>          1          2          3          4 
#> 100.000000  33.333333  11.111111   3.703704 
#> 

# Needs adjustment: alphas shrink at deeper levels
compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000,
                        max_depth = 5)
#>            1            2            3            4            5 
#> 0.0500000000 0.0125000000 0.0031250000 0.0007997540 0.0003946938 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>        1        2        3        4        5 
#>  1.00000  4.00000 15.62981 31.67012 20.97663 
#> 
#> attr(,"error_load")$sum_G
#> [1] 73.27656
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>         1         2         3         4         5 
#> 1.0000000 1.0000000 0.9768629 0.5065661 0.1655869 
#> 
#> attr(,"error_load")$critical_level
#> [1] 5
#> 
#> attr(,"error_load")$n_by_level
#>          1          2          3          4          5 
#> 1000.00000  250.00000   62.50000   15.62500    3.90625 
#> 
```
