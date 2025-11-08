# Compute local Simes p-value for a Vector of Child p-values

Given \\k\\ child p-values, computes the Simes p-value
\\\min\_{i=1\ldots k} \\ (k/i) \* p\_{(i)} \\\\, where \\p\_{(1)} \le
\ldots \le p\_{(k)}\\.

## Usage

``` r
local_simes(pvals_children, alpha = 0.05)
```

## Arguments

- pvals_children:

  Numeric vector of child p-values.

- alpha:

  Numeric scalar of alpha (not used in this function)

## Value

A single numeric value: the Simes combination p-value.

## Examples

``` r
local_simes(c(0.01, 0.04, 0.10, 0.20))
#> [1] 0.04
```
