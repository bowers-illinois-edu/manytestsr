# Handoff Summary

**Date:** 2026-03-01 **Branch:** `main` **Last commit:** `64af89e` —
Update HANDOFF.md for current session **Package version:** `0.0.4.1003`
(uncommitted — all changes pass
[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
with 0 errors, 0 warnings, 0 notes)

------------------------------------------------------------------------

## 1. Key Decisions Made

### Dependency lightening (REFACTOR_PLAN.md Steps 1–4)

- **8 packages moved from Imports to Suggests:** `stringi`, `tidygraph`,
  `ggraph`, `digest`, `ggplot2`, `Ckmeans.1d.dp`, `onlineFDR`, `hommel`.
  Each is now guarded with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and called
  via `pkg::fn()`.
- **`ClusterR` removed entirely.** The code that called `KMeans_rcpp`
  was already commented out; only the `@importFrom` tag remained.
- **`dataPreparation` replaced** with an internal `which_are_constant()`
  helper in `R/utils.R`. A one-liner:
  `which(vapply(df, function(x) length(unique(x)) < 2, logical(1)))`.
  Returns named integer indices. The only difference from
  [`dataPreparation::which_are_constant()`](https://rdrr.io/pkg/dataPreparation/man/which_are_constant.html)
  is that ours returns named indices (via
  [`which()`](https://rdrr.io/r/base/which.html)); tests use
  [`unname()`](https://rdrr.io/r/base/unname.html) for cross-validation.
- **Pipe and re-export breakage.** Removing `@import ggplot2` and
  `@import tidygraph` broke re-exports of `%>%`,
  [`filter()`](https://rdrr.io/r/stats/filter.html), `mutate()`, and
  [`grid::unit()`](https://rdrr.io/r/grid/unit.html). Fixed by replacing
  pipe chains with sequential explicit calls
  (`tidygraph::activate(g, "nodes"); g <- dplyr::mutate(g, ...)`) and
  prefixing [`grid::unit()`](https://rdrr.io/r/grid/unit.html).

### New functions: pPolyRank and pCombStephenson

- **[`pPolyRank()`](reference/pPolyRank.md)** — Tests Fisher’s sharp
  null using multiple polynomial rank scores simultaneously via
  [`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html).
  Computes within-block polynomial scores `(rank(Y) / (n_b + 1))^(r-1)`
  for each r in `r_vec` (default 2, 6, 10), passes them as a
  multivariate response. No new dependencies (uses existing `coin`
  import). Fast (~0.02s). Follows the exact same pattern as
  `pIndepDist`: builds `score1 + score2 + score3 ~ treatment | block`
  formula, uses `teststat = "quadratic"` by default.

- **[`pCombStephenson()`](reference/pCombStephenson.md)** — Wraps
  [`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html)
  for quantile-of-effects hypotheses. Default `k = n, c = 0` tests the
  sharp null (same as `pPolyRank` but via permutation-max combination,
  ~0.5s). For `k < n`, tests whether the k-th largest individual effect
  exceeds `c`. CMRSS is a Suggests dependency from GitHub.

### CMRSS k parameter semantics (critical for future work)

The `k` parameter in CMRSS indexes the k-th *largest* effect τ\_(k)
across all N units. Internally, the top `min(m, n-k)` treated units’
adjusted outcomes are set to infinity. **The test is degenerate when
`k ≤ n - m`** (all treated units get infinity → test statistic is
constant across permutations → p = 1 always). The wrapper warns on
degenerate k. At `k = n`, no units get infinity and the test reduces to
a standard stratified rank-sum.

### pPolyRank as the preferred sharp-null test

For the sharp null of no effects (which is what the tree-testing
procedure needs), `pPolyRank` is preferred over `pCombStephenson`
because: - 25x faster (asymptotic multivariate normal vs permutation
null) - No extra dependency (uses `coin`, already in Imports) - Produces
the same rank-sum statistics with polynomial scores - The “combined”
benefit (adaptive across r values) comes from coin’s multivariate test
infrastructure, which accounts for the correlation structure among score
functions via the joint covariance matrix

`pCombStephenson` remains valuable for quantile-of-effects hypotheses (k
\< n), where the CMRSS optimization over effect assignments is the core
contribution.

## 2. Files Changed and Why

### New files

| File                               | Purpose                                                                                       |
|------------------------------------|-----------------------------------------------------------------------------------------------|
| `R/utils.R`                        | Internal `which_are_constant()` replacing `dataPreparation`                                   |
| `tests/testthat/test_utils.R`      | 6 tests for `which_are_constant()`                                                            |
| `tests/testthat/test_stephenson.R` | 8 tests for [`pCombStephenson()`](reference/pCombStephenson.md) (skip if CMRSS not installed) |
| `tests/testthat/test_polyrank.R`   | 7 tests for [`pPolyRank()`](reference/pPolyRank.md)                                           |
| `man/pCombStephenson.Rd`           | Generated by roxygen2                                                                         |
| `man/pPolyRank.Rd`                 | Generated by roxygen2                                                                         |

### Modified files

| File                                       | What changed                                                                                                                                                                                                                                                                                                                                                                 |
|--------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `DESCRIPTION`                              | Version 0.0.4.1002 → 0.0.4.1003. Imports reduced to `Rcpp`, `Rfast`. Suggests gained `CMRSS`, `stringi`, `tidygraph`, `ggraph`, `digest`, `ggplot2`, `Ckmeans.1d.dp`, `onlineFDR`, `dataPreparation`, `hommel`. Remotes gained `bowers-illinois-edu/CMRSS`.                                                                                                                  |
| `NAMESPACE`                                | Regenerated. Reduced from ~125 to ~93 lines. Added `export(pCombStephenson)`, `export(pPolyRank)`. Removed all `importFrom` entries for moved packages.                                                                                                                                                                                                                      |
| `NEWS.md`                                  | Added 0.0.4.1003 entry documenting both new functions and dependency changes.                                                                                                                                                                                                                                                                                                |
| `.Rbuildignore`                            | Added `REFACTOR_PLAN\.md`                                                                                                                                                                                                                                                                                                                                                    |
| `R/pval_fns.R`                             | Removed 3 `@importFrom dataPreparation` tags, replaced 6 calls with internal helper. Added [`pPolyRank()`](reference/pPolyRank.md) (~100 lines) and [`pCombStephenson()`](reference/pCombStephenson.md) (~90 lines).                                                                                                                                                         |
| `R/test_statistics.R`                      | Replaced 1 [`dataPreparation::which_are_constant()`](https://rdrr.io/pkg/dataPreparation/man/which_are_constant.html) call.                                                                                                                                                                                                                                                  |
| `R/splitting_fns.R`                        | Removed `@importFrom ClusterR`, `@importFrom Ckmeans.1d.dp`, `@importFrom stringi`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards. Prefixed all calls with `pkg::`.                                                                                                                                                                               |
| `R/find_blocks.R`                          | Removed `@importFrom stringi`, `@importFrom digest`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards. Fixed dead `[KMeans_rcpp()]` doc link → `[Ckmeans.1d.dp::Ckmeans.1d.dp()]`.                                                                                                                                                                   |
| `R/reporting.R`                            | Removed `@import ggplot2`, `@import ggraph`, `@import tidygraph`, `@importFrom stringi`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards. Replaced all `%>%` pipes with sequential explicit calls. Replaced `unit()` with [`grid::unit()`](https://rdrr.io/r/grid/unit.html). Prefixed all calls with `pkg::`. Most complex changes in the session. |
| `R/alpha_fns.R`                            | Removed 3 `@importFrom onlineFDR`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards. Prefixed `Alpha_investing`, `SAFFRON`, `ADDIS` with `onlineFDR::`.                                                                                                                                                                                              |
| `R/p_adjustment.R`                         | Removed `@importFrom hommel`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guard. Prefixed `hommel()` with `hommel::`.                                                                                                                                                                                                                                  |
| `R/design_sensitivity.R`                   | Removed `@importFrom ggplot2`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guard. Prefixed all ggplot2 calls.                                                                                                                                                                                                                                          |
| `R/intersection_union_tests.R`             | Removed `@importFrom ggplot2`. Added [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guard. Prefixed all ggplot2 calls.                                                                                                                                                                                                                                          |
| `man/find_blocks.Rd`                       | Regenerated (doc link fix).                                                                                                                                                                                                                                                                                                                                                  |
| `tests/testthat/test_adaptations.R`        | Line 724: `stri_sub()` → [`substr()`](https://rdrr.io/r/base/substr.html) (stringi no longer imported in NAMESPACE).                                                                                                                                                                                                                                                         |
| `tests/testthat/test_reporting_plotting.R` | Lines 117–128, 171–182: Replaced bare `activate()` and `%>%` with [`tidygraph::activate()`](https://tidygraph.data-imaginist.com/reference/activate.html) explicit calls.                                                                                                                                                                                                    |

## 3. Current Blockers or Open Questions

- **No blockers.** All changes pass
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  with 0/0/0. Changes are uncommitted.
- **Open question: Should `pPolyRank` be the default `pfn` in
  [`find_blocks()`](reference/find_blocks.md)?** Currently
  [`find_blocks()`](reference/find_blocks.md) defaults to
  `pfn = pIndepDist`. The user may want to try `pfn = pPolyRank` for the
  tree search. This is a user decision, not a code issue.
- **Open question: Speed of `pCombStephenson` inside
  [`find_blocks()`](reference/find_blocks.md).** If `pCombStephenson` is
  used as `pfn` inside the tree search (hundreds of calls), the ~0.5s
  per call adds up. `pPolyRank` at ~0.02s per call is 25x faster for the
  sharp null. For quantile-of-effects hypotheses, there’s no fast
  alternative yet.

## 4. Important Context to Preserve

### CMRSS k parameter (critical footgun)

The `k` parameter in
[`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html)
indexes the k-th *largest* effect out of all N units. The test is
**degenerate** when `k ≤ n - m` because all treated units receive
`xi = Inf` internally, making the test statistic constant across
permutations (p = 1 always). The wrapper warns on this. At `k = n`, no
units get infinity and the test is a standard rank-sum.

### Pipe removal pattern

When `@import tidygraph` / `@import ggplot2` provided `%>%`, removing
those imports required replacing all pipe chains in `R/reporting.R` with
sequential explicit calls. The pattern used:

``` r
# Before:
res_graph %>% activate(nodes) %>% mutate(x = 1) %>% filter(y > 0)

# After:
res_graph <- tidygraph::activate(res_graph, "nodes")
res_graph <- dplyr::mutate(res_graph, x = 1)
res_graph <- dplyr::filter(res_graph, y > 0)
```

The same pattern was applied in
`tests/testthat/test_reporting_plotting.R`.

### Polynomial rank scores

The scores used by both `pPolyRank` and `pCombStephenson` are
`(rank(Y) / (n + 1))^(r-1)`: - r = 2: nearly linear in rank
(Wilcoxon-like, sensitive to location shifts) - r = 6, 10: increasingly
emphasize high ranks (sensitive to upper-tail effects)

### Test data outcomes

`make_test_data.R` creates three outcomes: - `Y`:
heterogeneous/canceling effects (positive in some blocks, negative in
others) - `Yhomog`: strong uniform positive effects across all blocks -
`Ynull`: no treatment effects

### Development workflow (CLAUDE_CODING.md)

Checkpoints: (1) write tests first, pause for review; (2) after
implementation, pause before running check; (3) after passing check,
pause for review.

## 5. What’s Done vs. What Remains

### Done

- All 11 steps of REFACTOR_PLAN.md implemented
- `which_are_constant()` internal helper created and tested
- 8 packages moved from Imports to Suggests with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guards
- `ClusterR` removed, `dataPreparation` replaced
- [`pCombStephenson()`](reference/pCombStephenson.md) wrapper created
  with correct `k = n` default for sharp null
- [`pPolyRank()`](reference/pPolyRank.md) created as fast coin-based
  alternative using polynomial rank scores
- Tests for all new functions (test_utils.R, test_stephenson.R,
  test_polyrank.R)
- All existing tests fixed for removed imports (pipes, activate,
  stri_sub, unit)
- Version bumped to 0.0.4.1003, NEWS.md updated
- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  passes with 0 errors, 0 warnings, 0 notes
- REFACTOR_PLAN.md added to .Rbuildignore

### Remains (not started)

- **Commit the changes.** All work is uncommitted. The user has not yet
  requested a commit.
- **Consider whether `pPolyRank` should be integrated into
  [`find_blocks()`](reference/find_blocks.md) as a `pfn` option** or
  documented as such in the vignette.
- **REFACTOR_PLAN.md cleanup.** The plan file is still in the repo root
  (ignored by .Rbuildignore). Could be deleted or moved once the user
  confirms the refactor is complete.
- **Potential future work:** A coin-based fast alternative for
  quantile-of-effects hypotheses (k \< n) if `pCombStephenson` speed
  becomes a bottleneck in tree search.
