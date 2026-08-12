# FIXME: stray object argument in the permutation branch of the p-value functions

**Status:** open. **Found:** 2026-06-22, while building Table 1 of the
NSF MMS proposal (`~/repos/power_for_policy/nsf_mms_2026`), which calls
`coin` directly partly to avoid this. **Severity:** real but narrow —
breaks only the `parallel = "snow"` backend of the Monte Carlo branch;
the default `parallel = "no"` and `parallel = "multicore"` work.
**Effort:** ~10 minutes plus a test.

This is a note for a future session working inside this package. Please
follow the package’s own rules in `CLAUDE.md`/`CLAUDE_CODING.md`: write
the test first, show it, then fix, then
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
(if any doc changes),
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html),
and bump the patch version in `DESCRIPTION`.

## What is wrong

Several p-value functions in `R/pval_fns.R` build the Monte Carlo
reference distribution like this (this is `pWilcox`, but the others are
identical):

``` r

    thedist <- coin::approximate(
      object,
      nresample = sims,
      parallel = parallel,
      ncpus = ncpu
    )
```

There is no variable named `object` in scope (confirmed: it is not in
the function arguments, not in the package namespace, not in `coin`). It
is a leftover. The call survives anyway because of how R matches
arguments. The signature is

    coin::approximate(nresample = 10000L, parallel = c("no","multicore","snow"),
                      ncpus = 1L, cl = NULL, B)

`nresample`, `parallel`, and `ncpus` are matched by name, so the
positional `object` falls through to the first remaining formal, which
is `cl` (the cluster object). Under `parallel = "no"` and
`parallel = "multicore"`, `cl` is never used, and R’s lazy evaluation
means the undefined `object` is never forced — so the call works and
returns a valid permutation p-value. Under `parallel = "snow"`, `coin`
uses `cl`, which forces `object`, and the call fails with:

    Error: object 'object' not found

So the permutation branch is fine in normal use (the default is
`parallel = "no"`) and silently fragile: it will break for any user who
asks for the snow backend, and the stray argument is confusing to read.

## Reproduce

``` r

library(manytestsr); data(example_dat)
sb <- subset(example_dat, blockF == "B080")   # n = 129
# simthresh above n forces the Monte Carlo branch:
pWilcox(sb, Y1 ~ trtF | blockF, simthresh = 200, sims = 200, parallel = "no")        # works
pWilcox(sb, Y1 ~ trtF | blockF, simthresh = 200, sims = 200, parallel = "multicore", ncpu = 2)  # works
pWilcox(sb, Y1 ~ trtF | blockF, simthresh = 200, sims = 200, parallel = "snow", ncpu = 2)        # ERROR: object 'object' not found
```

Note: the roxygen example in `pWilcox` passes `simthresh = 100`, but
block `B080` has 129 observations, so `nrow(dat) > simthresh` is TRUE
and the example takes the *asymptotic* branch — which is why
`example(pWilcox)` and `R CMD check` have not been catching this. Set
`simthresh` above the block size to exercise the Monte Carlo path.

## Where (all in `R/pval_fns.R`)

The identical `object,` line appears in five functions. Locations are
approximate (line numbers drift); search for the `coin::approximate(`
calls:

- `pOneway` (~ line 64)
- `pWilcox` (~ line 137)
- `pIndepDist` (~ line 291)
- `pTestTwice` (~ line 415)
- `pCombCauchyDist`(~ line 569)

`pPolyRank` (~ line 756) is **already correct** — its call has no
`object` and is the template for the fix:

``` r

    thedist <- coin::approximate(
      nresample = sims,
      parallel = parallel,
      ncpus = ncpu
    )
```

## The fix

Delete the `object,` line from each of the five calls so they match
`pPolyRank`. Nothing else changes; `nresample`/`parallel`/`ncpus` are
already named.

## Test to add first (tests/testthat)

The existing tests use `simthresh = 1` to force the *asymptotic* branch,
so the Monte Carlo branch is not covered. Add a test that exercises it
across all three backends for at least one affected function (and,
cheaply, for each of the five). For each: set `simthresh` above the
block size so the Monte Carlo branch runs, then assert the call returns
a single finite p-value in \[0, 1\] without error — in particular for
`parallel = "snow"`, which is the case that currently fails. A
regression test for `parallel = "snow"` is the one that matters; the
`"no"`/`"multicore"` cases guard against a fix that breaks them.

Keep the test fast (small `sims`, one small block). The statistical
principle being checked is that referring the statistic to its
randomization distribution returns a valid p-value regardless of the
(computational) parallel backend — the backend must not change validity.
