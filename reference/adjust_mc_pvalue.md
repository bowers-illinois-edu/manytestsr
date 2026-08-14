# Apply the Monte Carlo p-value convention

Converts the count/sims p-value that `coin` reports on the simulation
branch into \\(count + 1) / (sims + 1)\\. Asymptotic p-values are
returned unchanged.

## Usage

``` r
adjust_mc_pvalue(p, sims, thedist)
```

## Arguments

- p:

  Numeric p-value as returned by
  [`coin::pvalue()`](https://rdrr.io/pkg/coin/man/pvalue-methods.html).

- sims:

  Number of resamples passed to
  [`coin::approximate()`](https://rdrr.io/pkg/coin/man/NullDistribution.html).

- thedist:

  The distribution argument handed to the `coin` test: either the string
  `"asymptotic"` or an approximate null distribution. Used to decide
  whether the correction applies.

## Value

The corrected p-value, in \\\[1/(sims + 1), 1\]\\ on the simulation
branch and unchanged on the asymptotic branch.
