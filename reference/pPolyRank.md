# P-value function: Polynomial Rank Score Test

Tests Fisher's sharp null of no treatment effects using multiple
polynomial rank score functions simultaneously via
[`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html).
This provides adaptive sensitivity to effects at different parts of the
outcome distribution without pre-committing to a single scoring.

## Usage

``` r
pPolyRank(
  dat,
  fmla = Y ~ trtF | blockF,
  r_vec = c(2, 6, 10),
  teststat = "quadratic",
  simthresh = 20,
  sims = 1000,
  parallel = "no",
  ncpu = NULL
)
```

## Arguments

- dat:

  An object inheriting from class data.frame (will be coerced to
  data.table internally if needed).

- fmla:

  A formula: outcome ~ treatment \| block.

- r_vec:

  Numeric vector. Polynomial rank score parameters. Default
  `c(2, 6, 10)`. r = 2 gives Wilcoxon-like scores; larger r emphasizes
  extreme ranks.

- teststat:

  Character. Type of test statistic passed to
  [`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html).
  `"quadratic"` (default) gives an omnibus chi-squared statistic that
  pools across score functions; `"maximum"` takes the max of the
  standardized statistics across score functions.

- simthresh:

  Integer. Below this number of observations, use permutation-based
  inference instead of asymptotic approximation.

- sims:

  Integer. Number of permutations when using simulation-based inference.

- parallel:

  Character. `"no"` (default), `"snow"`, or `"multicore"`, passed to
  [`coin::approximate()`](https://rdrr.io/pkg/coin/man/NullDistribution.html).

- ncpu:

  Integer. Number of CPUs for parallel operation.

## Value

A p-value (numeric scalar).

## Details

For each value of r in `r_vec`, the function computes within-block
polynomial rank scores: `(rank(Y) / (n_b + 1))^(r-1)` where `n_b` is the
block size. These scores are passed as a multivariate response to
[`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html),
which computes a joint test statistic that accounts for the correlation
structure among the score functions.

With r = 2 the scores are nearly linear in rank (Wilcoxon-like). Larger
r values place increasing weight on observations with high ranks,
providing sensitivity to treatment effects concentrated in the upper
tail of the outcome distribution.

This function tests the same sharp null as
[`pOneway`](https://bowers-illinois-edu.github.io/manytestsr/reference/pOneway.md)
and
[`pIndepDist`](https://bowers-illinois-edu.github.io/manytestsr/reference/pIndepDist.md)
but uses polynomial rank scores from the combined Stephenson rank test
framework (Kim, Li, and Bowers). For quantile-of-effects hypotheses
(whether the k-th largest individual effect exceeds a threshold), see
[`pCombStephenson`](https://bowers-illinois-edu.github.io/manytestsr/reference/pCombStephenson.md).

## Examples

``` r
# Example using built-in data
data(example_dat, package = "manytestsr")
library(data.table)

# Test for treatment effect using polynomial rank scores
single_block <- as.data.table(subset(example_dat, blockF == "B080"))
p_val <- pPolyRank(single_block, Y1 ~ trtF | blockF)
print(p_val)
#> [1] 0.627839
```
