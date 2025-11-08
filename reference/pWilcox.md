# P-value function: Wilcox Test

These functions accept a data frame and perhaps test specific arguments
(like whether or not the test will be asympotic or simulation based). It
produces a p-value.

## Usage

``` r
pWilcox(
  dat,
  fmla = YContNorm ~ trtF | blockF,
  simthresh = 20,
  sims = 1000,
  parallel = "no",
  ncpu = NULL
)
```

## Arguments

- dat:

  An object inheriting from class data.frame

- fmla:

  outcome~treatment factor \| block factor (following `coin` API).

- simthresh:

  Below which number of total observations should the p-value functions
  use permutations rather than asymptotic approximations

- sims:

  Either NULL (meaning use an asymptotic reference dist) or a number
  (meaning sampling from the randomization distribution implied by the
  formula)

- parallel:

  Should the function use multicore processing for permutation based
  testing. Default is no. But could be "snow" or "multicore" following
  `approximate` in the coin package.

- ncpu:

  is the number of workers (for "snow") or cores (for "multicore").

## Value

A p-value

## Examples

``` r
# Example using Wilcoxon rank-sum test
data(example_dat, package = "manytestsr")

# Test for treatment effect on Y1 within a single block
single_block <- subset(example_dat, blockF == "B080")
p_val <- pWilcox(single_block, Y1 ~ trtF | blockF, parallel = "no")
print(p_val)
#> [1] 0.2547046

# Compare with permutation-based version
p_val_perm <- pWilcox(single_block, Y1 ~ trtF | blockF,
  simthresh = 100, sims = 500, parallel = "no"
)
print(p_val_perm)
#> [1] 0.2547046
```
