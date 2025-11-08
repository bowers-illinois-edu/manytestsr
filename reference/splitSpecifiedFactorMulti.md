# A set of pre-specified splits

This function allows for more than two splits at each level

## Usage

``` r
splitSpecifiedFactorMulti(bid, x)
```

## Arguments

- bid:

  Block id

- x:

  Is a a factor with levels like "state.district.school". The splits
  will occur from left to right depending on whether there is existing
  variation at that level
