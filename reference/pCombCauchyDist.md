# P-value function: Cauchy Combined Indepence Test

P-value function: Cauchy Combined Indepence Test

## Usage

``` r
pCombCauchyDist(
  dat,
  fmla = YcontNorm ~ trtF | blockF,
  simthresh = 20,
  sims = 1000,
  parallel = "no",
  ncpu = NULL,
  distfn = fast_dists_and_trans_nomax_hybrid
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
  matrix) of the same number of rows as the dat. Here until further
  notice users should leave it at `mean_dist_raw_rank_tanh` since that
  function is purpose built for this function.

## Value

A p-value

## Details

This function combines the p-values from univariate tests using using
(1) the raw outcome, (2) the rank-transformed outcome, (3) the tanh
transformed raw outcome (another test statistic that is more sensitive
when there is skew), (4) the mean difference in raw outcome Euclidean
distances bwtween treated and control observations within block; (5) the
mean difference in ranked outcome Euclidean distances bwtween treated
and control observations; (6) a Hotelling-T style quadratic combination
of the preceding test statistics. Inspired by Rizzo and Székely's work
on Euclidean distance based testing and by Liu and Xie (2020) on the
Cauchy Combination Test and Hansen and Bowers (2008) on omnibus tests.
Distance and ranks and other transformations are all calculated by block
when block is supplied in the formula.
