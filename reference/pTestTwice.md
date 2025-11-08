# P-value function: Testing twice

These functions accept a data frame and perhaps test specific arguments
(like whether or not the test will be asymptotic or simulation based).
It produces a p-value.

## Usage

``` r
pTestTwice(
  dat,
  fmla = YcontNorm ~ trtF | blockF,
  simthresh = 20,
  sims = 1000,
  parallel = "yes",
  ncpu = NULL,
  groups = NULL
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

- groups:

  Currently unused parameter, reserved for future functionality

## Value

A p-value

## Details

For now, this function does an omnibus-style max-T test using (1) the
raw outcome and (2) a rank transformed raw outcome. Inspired by
Rosenbaum (2008) on Testing Twicee
