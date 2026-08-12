# Compute Adaptive Alpha Levels by Tree Depth

Computes the per-depth significance levels for a regular k-ary tree
under the adaptive-alpha framework of Appendix B of the supplement.
First checks whether natural gating suffices (total error load \\\le
1\\); if so, returns nominal alpha at every level. Otherwise, computes
adjusted significance levels that compensate for the error load at each
depth, optionally with user-supplied depth-wise budget weights.

## Usage

``` r
compute_adaptive_alphas(
  k,
  delta_hat,
  N_total,
  max_depth = 20L,
  thealpha = 0.05,
  budget_weights = NULL
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

- budget_weights:

  Controls how the error budget is allocated across depths. Same options
  as in
  [`compute_adaptive_alphas_tree`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas_tree.md):
  `NULL` (default, telescoping), `"equal"`, `"proportional"`, or a
  numeric vector of length `max_depth - 1` that sums to at most 1.

## Value

Named numeric vector of adjusted alpha levels, one per depth (1 through
`max_depth`). Names are depth levels as characters. Has attribute
`"error_load"` containing the
[`compute_error_load`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_error_load.md)
result, so the caller can inspect whether adjustment was needed.

## Details

The function first calls
[`compute_error_load`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_error_load.md)
to assess whether natural gating suffices. When \\\sum G\_\ell \le 1\\,
no adjustment is needed and nominal `thealpha` is returned at every
level (and `budget_weights` is ignored, since the tree protects itself).

When adjustment is needed and `budget_weights = NULL` (per-depth
telescoping), the formula at level \\\ell\\ is: \$\$\alpha\_\ell^{adj} =
\min\left\\\alpha,\\ \frac{\alpha}{k^{(\ell-1)} \cdot
\prod\_{j=1}^{\ell-1} \hat\theta_j} \right\\\$\$

When `budget_weights` is supplied, the formula becomes
\\\alpha\_\ell^{adj} = \min\\\alpha, w\_\ell \alpha / G\_\ell\\\\, where
\\w\_\ell\\ is the resolved weight at depth \\\ell\\. The constraint
\\\sum w\_\ell \le 1\\ is what gives the FWER guarantee via the
budget-weighted theorem in the supplement.

The FWER guarantee requires that power is not underestimated (i.e.,
\\\hat\theta_j \geq \theta_j\\). In practice this means using a
conservatively large `delta_hat`.

## Examples

``` r
# Natural gating sufficient: all alphas = 0.05
compute_adaptive_alphas(k = 3, delta_hat = 0.2, N_total = 100,
                        max_depth = 4)
#>          1          2          3          4 
#> 0.05000000 0.03229941 0.05000000 0.05000000 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>         1         2         3         4 
#> 0.0000000 1.5480158 0.9810764 0.3009919 
#> 
#> attr(,"error_load")$sum_G
#> [1] 2.830084
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>          1          2          3          4 
#> 0.51600527 0.21125461 0.10226587 0.06713787 
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
#> 0.0500000000 0.0125000000 0.0031250000 0.0007997540 0.0003946616 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>         1         2         3         4         5 
#>   0.00000   4.00000  16.00000  62.51922 126.69082 
#> 
#> attr(,"error_load")$sum_G
#> [1] 209.21
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>         1         2         3         4         5 
#> 1.0000000 1.0000000 0.9768629 0.5066075 0.1671852 
#> 
#> attr(,"error_load")$critical_level
#> [1] 5
#> 
#> attr(,"error_load")$n_by_level
#>          1          2          3          4          5 
#> 1000.00000  250.00000   62.50000   15.62500    3.90625 
#> 

# Budget-weighted: equal allocation across depths 2..max_depth
compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000,
                        max_depth = 5, budget_weights = "equal")
#>            1            2            3            4            5 
#> 5.000000e-02 3.125000e-03 7.812500e-04 1.999385e-04 9.866539e-05 
#> attr(,"error_load")
#> attr(,"error_load")$G
#>         1         2         3         4         5 
#>   0.00000   4.00000  16.00000  62.51922 126.69082 
#> 
#> attr(,"error_load")$sum_G
#> [1] 209.21
#> 
#> attr(,"error_load")$needs_adjustment
#> [1] TRUE
#> 
#> attr(,"error_load")$thetas
#>         1         2         3         4         5 
#> 1.0000000 1.0000000 0.9768629 0.5066075 0.1671852 
#> 
#> attr(,"error_load")$critical_level
#> [1] 5
#> 
#> attr(,"error_load")$n_by_level
#>          1          2          3          4          5 
#> 1000.00000  250.00000   62.50000   15.62500    3.90625 
#> 
```
