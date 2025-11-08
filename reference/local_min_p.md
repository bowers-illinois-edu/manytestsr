# Unadjusted local minimal p-value

Given \\k\\ child p-values, returns the highest p-value below alpha; if
none are below alpha returns the smallest p-value.

## Usage

``` r
local_min_p(pvals_children, alpha = 0.05)
```

## Arguments

- pvals_children:

  Numeric vector of child p-values.

- alpha:

  Numeric scalar of alpha

## Value

A single numeric value.

## Examples

``` r
local_min_p(c(0.01, 0.04, 0.10, 0.20))
#> [1] 0.04
local_min_p(c(0.10, 0.20))
#> [1] 0.1
```
