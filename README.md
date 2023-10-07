# A package to do many tests to localize causal effects in (sets of) experimental blocks

## Development Info

```
devtools::document()
devtools::check()
```

## Notes:

We are using three main C functions: `fast_dists_by_unit_arma2_par` if we are
using openmp (`parallel="yes"`), `fast_dists_and_trans` for small n, and
`fast_dists_and_trans_by_unit_arma` for N>20. (Where 20 is chosen without a lot
of profiling. Basic idea is not to hold large matrices in memory. Matrices are
faster when N is not to large, but bog down when we have large number of
observations.
