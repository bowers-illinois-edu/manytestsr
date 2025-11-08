# Compute local Hommel p-value for a Vector of Child p-values

Given \\k\\ child p-values, computes the Hommel-adjusted p-values

## Usage

``` r
local_hommel_all_ps(pvals_children, alpha = 0.05)
```

## Arguments

- pvals_children:

  Numeric vector of child p-values.

- alpha:

  Numeric scalar of alpha (not used in this function)

## Value

A vector of adjusted p-values

## Examples

``` r
local_hommel_all_ps(c(0.01, 0.04, 0.10, 0.20))
#> [1] 0.04 0.12 0.20 0.20
```
