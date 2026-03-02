# P-value function: Combined Stephenson Rank Test

Wraps
[`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html)
to provide a consistent formula interface for the combined Stephenson
rank test of Kim, Li, and Bowers. Tests quantile-of-effects hypotheses:
whether the k-th largest individual treatment effect exceeds a threshold
c.

## Usage

``` r
pCombStephenson(
  dat,
  fmla = Y ~ trtF | blockF,
  k = NULL,
  c = 0,
  r_vec = c(2, 6, 10),
  weight_name = "asymp.opt",
  null_max = 10^5,
  opt_method = "ILP_auto",
  comb_method = 1
)
```

## Arguments

- dat:

  An object inheriting from class data.frame (will be coerced to
  data.table internally if needed).

- fmla:

  A formula: outcome ~ treatment \| block.

- k:

  Integer. The quantile index: test whether the k-th largest individual
  effect exceeds `c`. Default is `n` (the total number of units), which
  tests Fisher's sharp null of no effects. Must satisfy `k > n - m` for
  the test to be non-trivial. Smaller k tests quantile-of-effects
  hypotheses: whether at least `n - k + 1` units have effects exceeding
  `c`.

- c:

  Numeric. The threshold for the null hypothesis. Default 0.

- r_vec:

  Integer vector. Polynomial rank score tuning parameters. Default
  `c(2, 6, 10)`. These control the weight placed on different parts of
  the rank distribution: r = 2 is close to Wilcoxon, larger r emphasizes
  extreme ranks.

- weight_name:

  Character. Weighting scheme across blocks. `"asymp.opt"` (default)
  uses asymptotically optimal weights; `"dis.free"` uses
  distribution-free weights.

- null_max:

  Integer. Number of permutations for the randomization null
  distribution. Default 100000.

- opt_method:

  Character. Optimization method for the test statistic. Default
  `"ILP_auto"` (integer linear programming with automatic solver
  selection). See
  [`?CMRSS::pval_comb_block`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html)
  for options.

- comb_method:

  Integer. 1 = aggregate across strata then combine across methods
  (default); 2 = combine across methods within each stratum then
  aggregate (often more powerful).

## Value

A p-value (numeric scalar).

## Details

The combined test uses polynomial rank scores at multiple tuning
parameters (controlled by `r_vec`) and takes the maximum of the
standardized statistics, avoiding the need to choose a single tuning
parameter. This wrapper builds the `methods.list.all` configuration that
[`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html)
expects.

The `k` parameter indexes the k-th largest treatment effect, where
tau\_(1) \>= tau\_(2) \>= ... \>= tau\_(n). The test asks whether
tau\_(k) exceeds `c`. Internally, CMRSS sets the top `min(m, n-k)`
treated units' adjusted outcomes to infinity. The test is non-trivial
only when `k > n - m` (where n is total units and m is number treated),
because otherwise all treated units receive infinity and the test
statistic is constant across all permutations.

The default `k = n` and `c = 0` tests Fisher's sharp null of no effects
for any unit. At k = n no treated units receive infinity, so the test
statistic reduces to a standard stratified rank-sum statistic — but
combined across multiple polynomial score functions (controlled by
`r_vec`), which provides adaptive sensitivity without pre-committing to
a single scoring. Smaller k values test quantile-of-effects hypotheses:
whether at least `n - k + 1` units have effects exceeding `c`.

## Examples

``` r
# \donttest{
# Requires CMRSS package:
# remotes::install_github("bowers-illinois-edu/CMRSS")
if (requireNamespace("CMRSS", quietly = TRUE)) {
  data(example_dat, package = "manytestsr")
  library(data.table)
  idat <- as.data.table(example_dat)
  p <- pCombStephenson(idat, Y1 ~ trtF | blockF)
  print(p)
}
#> Error in get_default_solver(): No solver available. Please install either 'highs' (recommended, open-source) or 'gurobi'.
#> To install highs: install.packages('highs')
#> For gurobi, see: https://www.gurobi.com/documentation/current/quickstart_mac/r_ins_the_r_package.html
# }
```
