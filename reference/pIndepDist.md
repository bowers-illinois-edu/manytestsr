# P-value function: Independence Treatment Distance Test

These functions accept a data frame and perhaps test specific arguments
(like whether or not the test will be asymptotic or simulation based).
It produces a p-value.

## Usage

``` r
pIndepDist(
  dat,
  fmla = YcontNorm ~ trtF | blockF,
  simthresh = 20,
  sims = 1000,
  parallel = "yes",
  ncpu = NULL,
  distfn = fast_dists_and_trans_hybrid
)
```

## Arguments

- dat:

  An object inheriting from class data.frame

- fmla:

  A formula appropriate to the function. Here it should be something
  like outcome~treatment\|block

- simthresh:

  is the size of the data below which we use direct permutations for
  p-values

- sims:

  Either NULL (meaning use an asymptotic reference dist) or a number
  (meaning sampling from the randomization distribution implied by the
  formula)

- parallel:

  is "no" then parallelization is not required, otherwise it is
  "multicore" or "snow" in the call to
  [`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html)
  (see help for coin::approximate()). Also, if parallel is not "no" and
  `adaptive_dist_function` is TRUE, then an openmp version of the
  distance creation function is called using `ncpu` threads (or
  `parallel::detectCores(logical=FALSE)` cores).

- ncpu:

  is number of cpus to be used for parallel operation.

- distfn:

  is a function that produces one or more vectors (a data frame or
  matrix) of the same number of rows as the dat

## Value

A p-value

## Details

For now, this function does an omnibus-style chi-square test using (1)
the ratio of distances to controls to distances to treated observations
within block; (2) the rank of distances to controls for each unit; and
(3) the raw outcome.

Although the distances are calculated by block, our profiling suggests
that it is better to parallelize the distance creation `distfn` (done
here in C++ in the `fastfns.cpp` file) rather than use the `data.table`
approach of
[`setDTthreads()`](https://rdrr.io/pkg/data.table/man/openmp-utils.html).
So, here we assume that the threads for data.table are 1.

## Examples

``` r
# \donttest{
# Example using distance-based independence test
data(example_dat, package = "manytestsr")
library(data.table)

# Test for treatment effect using distance-based approach
single_block <- as.data.table(subset(example_dat, blockF == "B080"))
p_val <- pIndepDist(single_block, Y1 ~ trtF | blockF, parallel = "no")
print(p_val)
#> [1] 0.5394499

# Test with different outcome variable
p_val2 <- pIndepDist(single_block, Y2 ~ trtF | blockF, parallel = "no")
print(p_val2)
#> [1] 0.7764711
# }
```
