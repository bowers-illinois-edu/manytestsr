# Unadjusted local step (pass-through)

Returns the input p-values unmodified.

## Usage

``` r
local_unadj_all_ps(pvals_children, alpha = 0.05)
```

## Arguments

- pvals_children:

  Numeric vector of child p-values.

- alpha:

  Numeric scalar of alpha (not used)

## Value

A numeric vector: the same as input.

## Examples

``` r
local_unadj_all_ps(c(0.01, 0.04))
#> [1] 0.01 0.04
```
