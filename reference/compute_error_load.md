# Compute Error Load for Natural Gating Assessment

Computes the error load at each tree level, which measures the expected
number of all-null sibling groups tested. When the total error load is
at most 1, the unadjusted procedure (fixed alpha at every node) controls
FWER at level alpha — no adjustment needed. When it exceeds 1, adaptive
alpha adjustment is required.

## Usage

``` r
compute_error_load(
  k = NULL,
  delta_hat,
  N_total = NULL,
  node_dat = NULL,
  max_depth = 20L,
  thealpha = 0.05
)
```

## Arguments

- k:

  Branching factor. Either a scalar (constant k at all levels) or an
  integer vector of length `max_depth` where `k[ell]` is the branching
  factor at level `ell`. Used in parametric mode; ignored when
  `node_dat` is provided.

- delta_hat:

  Estimated standardized effect size (e.g., Cohen's d). Used to compute
  power at each level via the normal approximation.

- N_total:

  Total sample size at the root level. Used in parametric mode; ignored
  when `node_dat` is provided.

- node_dat:

  Optional data.frame or data.table with columns `nodenum`, `parent`,
  `depth`, and `nodesize`. When provided, the function computes per-node
  power from `nodesize` and aggregates error load by depth. This
  supports irregular trees (e.g., DPP design with unequal splits).

- max_depth:

  Maximum tree depth to compute. In parametric mode, defaults to 20. In
  tree mode, inferred from `node_dat`.

- thealpha:

  Nominal significance level (default 0.05).

## Value

A named list with components:

- `G`:

  Numeric vector of error load at each level (length `max_depth`).

- `sum_G`:

  Total error load: `sum(G)`.

- `needs_adjustment`:

  Logical: `TRUE` when `sum_G > 1`.

- `thetas`:

  Estimated power (conditional rejection probability) at each level.

- `critical_level`:

  First level where theta \< 1/k (where natural gating begins to
  dominate), or `NA` if theta never drops below 1/k.

- `n_by_level`:

  Sample size at each level.

## Details

Two interfaces are provided. The **parametric** interface (`k`,
`delta_hat`, `N_total`) assumes a regular k-ary tree with equal
splitting: sample size at level ell is \\N / k^{ell-1}\\. The **tree**
interface (`node_dat`) accepts per-node sample sizes from an actual
(possibly irregular) tree, as returned by
[`find_blocks`](find_blocks.md).

In **parametric mode**, the error load at level ell is: \$\$G\_\ell =
k^{\ell-1} \prod\_{j=0}^{\ell-1} \theta_j\$\$ where \\\theta_j =
\Phi(\delta \sqrt{n_j} - z\_{1-\alpha/2})\\ and \\n_j = N / k^{j}\\.

In **tree mode**, the function computes per-node power from the actual
sample sizes. For each node at depth d, it computes \\\theta =
\Phi(\delta \sqrt{n\_{\text{node}}} - z\_{1-\alpha/2})\\. The error load
at depth d is the sum over all nodes at depth d of the product of theta
values along the path from root to that node's parent, which generalizes
the regular-tree formula to irregular branching.

## Examples

``` r
# Parametric mode: regular 3-ary tree, moderate effect
compute_error_load(k = 3, delta_hat = 0.3, N_total = 500)
#> $G
#>            1            2            3            4            5            6 
#> 9.999990e-01 2.916379e+00 5.326176e+00 4.022713e+00 1.354743e+00 2.562653e-01 
#>            7            8            9           10           11           12 
#> 3.343749e-02 3.475357e-03 3.153782e-04 2.642414e-05 2.113216e-06 1.644896e-07 
#>           13           14           15           16           17           18 
#> 1.260452e-08 9.571413e-10 7.230190e-11 5.445132e-12 4.093620e-13 3.074454e-14 
#>           19           20 
#> 2.307679e-15 1.731555e-16 
#> 
#> $sum_G
#> [1] 14.91353
#> 
#> $needs_adjustment
#> [1] TRUE
#> 
#> $thetas
#>          1          2          3          4          5          6          7 
#> 0.99999897 0.97212721 0.60876590 0.25175746 0.11225782 0.06305386 0.04349332 
#>          8          9         10         11         12         13         14 
#> 0.03464531 0.03024900 0.02792852 0.02665765 0.02594617 0.02554268 0.02531212 
#>         15         16         17         18         19         20 
#> 0.02517981 0.02510368 0.02505982 0.02503452 0.02501993 0.02501150 
#> 
#> $critical_level
#> [1] 4
#> 
#> $n_by_level
#>            1            2            3            4            5            6 
#> 5.000000e+02 1.666667e+02 5.555556e+01 1.851852e+01 6.172840e+00 2.057613e+00 
#>            7            8            9           10           11           12 
#> 6.858711e-01 2.286237e-01 7.620790e-02 2.540263e-02 8.467544e-03 2.822515e-03 
#>           13           14           15           16           17           18 
#> 9.408382e-04 3.136127e-04 1.045376e-04 3.484586e-05 1.161529e-05 3.871762e-06 
#>           19           20 
#> 1.290587e-06 4.301958e-07 
#> 

# High power, wide tree: needs adjustment
res <- compute_error_load(k = 10, delta_hat = 0.5, N_total = 5000,
                          max_depth = 3)
res$sum_G     # likely > 1
#> [1] 105.2438
res$needs_adjustment
#> [1] TRUE
```
