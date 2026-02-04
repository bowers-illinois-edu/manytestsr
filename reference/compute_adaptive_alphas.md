# Compute Adaptive Alpha Levels by Tree Depth

Implements Algorithm 1 from the adaptive alpha appendix. Computes
adjusted significance levels for each tree level based on estimated
power decay through the tree.

## Usage

``` r
compute_adaptive_alphas(
  k,
  delta_hat,
  N_total,
  tau = 0.1,
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

- tau:

  Cumulative power threshold (default 0.1). When cumulative power drops
  below `tau`, natural gating is deemed sufficient and nominal alpha is
  used.

- max_depth:

  Maximum depth to compute (default 20).

- thealpha:

  Nominal significance level (default 0.05).

## Value

Named numeric vector of adjusted alpha levels, one per depth (1 through
`max_depth`). Names are depth levels as characters.

## Details

The formula at level \\\ell\\ is: \$\$\alpha\_\ell^{adj} = \alpha /
(k^{(\ell-1)} \cdot \prod\_{j=1}^{\ell-1} \hat\theta_j)\$\$ when
cumulative power exceeds `tau`, and \\\alpha\_\ell^{adj} = \alpha\\
otherwise.

The FWER guarantee (Theorem in Appendix D) requires that power estimates
are not overestimated. In practice this means using a conservatively
large `delta_hat`.

## Examples

``` r
# Alpha schedule for a 4-ary tree with moderate effect
compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 1000)
#>           1           2           3           4           5           6 
#> 0.050000000 0.012500000 0.003125000 0.000799754 0.050000000 0.050000000 
#>           7           8           9          10          11          12 
#> 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 
#>          13          14          15          16          17          18 
#> 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 0.050000000 
#>          19          20 
#> 0.050000000 0.050000000 

# Binary tree with very high power — approaches Bonferroni
compute_adaptive_alphas(k = 2, delta_hat = 0.5, N_total = 1e6,
                        tau = 0, max_depth = 5)
#>        1        2        3        4        5 
#> 0.050000 0.025000 0.012500 0.006250 0.003125 
```
