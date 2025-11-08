# BH local adjustment

Performs Benjamini-Hochberg adjustment on a vector of p-values.

## Usage

``` r
local_bh_all_ps(pvals_children, alpha = 0.05)
```

## Arguments

- pvals_children:

  Numeric vector of child p-values.

- alpha:

  Numeric scalar of alpha (not used)

## Value

A numeric vector of BH-adjusted p-values.

## Examples

``` r
local_bh_all_ps(c(0.01, 0.04, 0.10, 0.20))
#> [1] 0.0400000 0.0800000 0.1333333 0.2000000
```
