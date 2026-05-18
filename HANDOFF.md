# Handoff Summary

**Date:** 2026-05-18
**Branch:** `main` (1 commit ahead of `origin/main`, not pushed)
**Last commit:** `cb3604e` --- Fix pCombStephenson for CMRSS 0.2.6+ and switch deps to pak
**Package version:** `0.0.4.1005`
**Working tree:** Pre-existing alpha-adaptive changes still uncommitted (see section 2)
**`make check` status:** **0 errors, 0 warnings, 0 notes** (clean, ~1m 37s, vignettes built)

---

## 1. Key Decisions Made

### Switched Makefile dependency install from devtools to pak

User asked whether pak is the new recommended install approach. Confirmed:
pak is developed at Posit by Gabor Csardi --- complement, not competitor, to
devtools. *R Packages* 2nd ed. teaches pak first for installation. pak's
advantages here: parallel downloads, native handling of `Remotes:` (the three
GitHub deps in this package), better system-requirement messaging, faster
solver. Fully compatible with Makefile workflows --- pak is just a function
call replacing `devtools::install_deps`.

User chose: pak for the `dependencies` target; devtools stays for `check`,
`test`, `document`, `build`.

### Made `pCombStephenson` compatible with CMRSS 0.2.6+

Diagnosed 6 test failures in `test_stephenson.R` as a CMRSS upstream API
change. Confirmed by reading `CMRSS::pval_comb_block` source. Two distinct
changes upstream:

1. **CMRSS internal `k` is now constrained to 1..m (not 1..n).** Source:
   `p <- m - k; if (k < 1 || k > m) stop(...)`. The internal LP coverage
   constraint is `p = m - k` where `p` is the number of treated units whose
   adjusted outcomes get set to Inf.
2. **CMRSS now always calls `solve_optimization()`** even for the sharp-null
   path. Requires `highs` (open-source) or `gurobi`.

Semantic translation table:

| Concept                          | Old CMRSS (paper notation)    | New CMRSS 0.2.6+            |
|----------------------------------|-------------------------------|-----------------------------|
| Rank index range                 | k in 1..n                     | k in 1..m                   |
| Number of treated set to Inf     | min(m, n - k)                 | m - k                       |
| Sharp null (no Inf)              | k = n                         | k = m                       |
| Degenerate (all Inf, p = 1)      | k <= n - m                    | unreachable (k >= 1)        |
| Mapping new <- old               | ---                           | `cmrss_k = k - (n - m)`     |

User chose: keep paper-notation `k` in 1..n at the wrapper's user-facing API
(preserves existing call sites and paper conventions), translate internally
to `cmrss_k = k - (n - m)` before calling CMRSS. The degeneracy path
(`k <= n - m`) now warns and returns `p = 1` *without* calling CMRSS, avoiding
the LP-solver requirement on that branch. Installed `highs` and added to
Suggests.

### Vignettes on by default

Quarto CLI is installed locally (`/usr/local/bin/quarto`, v1.9.37), so the
default `make check` builds vignettes. Added `make check-fast` as an escape
hatch (skips vignettes and manual). Added `make check-cran` for stricter
pre-submission checking.

### Versioning and commit scope

This session's work claimed the existing `0.0.4.1005` version bump (already
in the working tree from the prior uncommitted session) and added a NEWS
entry for 1005 covering the CMRSS-compat fix and `highs` dependency. The
alpha-adaptive uncommitted work, when it eventually gets committed, will
need to bump to `0.0.4.1006` and add its own NEWS section.

Commit `cb3604e` staged only this session's files. The pre-existing
alpha-adaptive working-tree changes were left untouched.

## 2. Files Changed and Why

### Committed this session (`cb3604e`)

| File | What changed |
|------|-------------|
| `Makefile` | (1) Rewrote `dependencies` target to use `pak::pkg_install('.', dependencies = TRUE, upgrade = FALSE, ask = FALSE)`. (2) Added `check-fast` target (`vignettes = FALSE, manual = FALSE, --no-build-vignettes --no-manual`). (3) Added `check-cran` target (`cran = TRUE, remote = TRUE, manual = TRUE`). (4) Added header comment noting first-time setup needs `make dependencies`. Other targets (`interactive`, `test`, `check`, `document`, `build`, `clean`, `spell-check-DESCRIPTION`) unchanged. |
| `R/pval_fns.R` | Patched `pCombStephenson()`. (a) Added `if (k > n) stop(...)` validation. (b) Degenerate path (`k <= n - m`) now warns *and returns `p = 1`* without calling CMRSS. (c) Added the translation `cmrss_k <- k - (n - m)` before the CMRSS call, and passes `cmrss_k` instead of `k`. (d) Updated `@details` text to describe the new internal mapping and the `highs`/`gurobi` LP-solver requirement. Default `k = n` (sharp null) unchanged at the user-facing API. |
| `tests/testthat/test_stephenson.R` | Added `skip_if_not_installed("highs")` to the 4 tests that exercise the LP path. The "warns on degenerate k" test is unchanged because the wrapper now short-circuits before CMRSS. Updated the "matches direct CMRSS call" test to compute `cmrss_k = m` instead of `k = n` for the direct call. |
| `DESCRIPTION` | Added `highs` to `Suggests:` between `CMRSS` and `stringi`. Version `0.0.4.1005` was already bumped from `1004` in the pre-existing working tree; this commit claims that bump for the CMRSS-compat work. |
| `NEWS.md` | Added `# manytestsr 0.0.4.1005` section above the existing `1003` section. Covers the CMRSS 0.2.6+ compatibility bug fix and the new `highs` Suggests. |
| `HANDOFF.md` | This file. |
| `MEMORY.md` (auto-memory, not in repo) | Updated `pCombStephenson` description with new CMRSS API and `highs` requirement. |

### Still uncommitted (pre-existing, NOT touched this session)

These were modified/untracked before this session began (since the
2026-04-08 handoff, which incorrectly claimed the tree was clean). Their
content and intent are not documented anywhere. **The next session should
read these diffs and the two new test files to understand intent before
committing.**

```
M R/alpha_adaptive.R
M man/alpha_adaptive.Rd
M man/alpha_adaptive_tree_pruned.Rd
M man/compute_adaptive_alphas.Rd
M tests/testthat/test_alpha_adaptive.R
M tests/testthat/test_alpha_adaptive_pruned.R
M tests/testthat/test_alpha_adaptive_tree.R
?? tests/testthat/test_irregular_tree_warning.R
?? tests/testthat/test_weak_fwer_global_null.R
```

All of these pass `make check` cleanly (verified after `cb3604e`), so the
prior work is at least internally consistent. But the intent is undocumented
and they have not been committed.

### Recent commits

| Commit | What |
|--------|------|
| `cb3604e` | **This session.** Fix `pCombStephenson` for CMRSS 0.2.6+ and switch deps to pak. |
| `3463e0c` | Added `test_budget_weights.R` (18 tests) and `test_switching_corollary.R` (5 tests) |
| `1abff53` | Extended `R/alpha_adaptive.R` and `R/pval_fns.R` with budget-weight and switching corollary support |
| `9abeb8e` | Dependency lightening, `pPolyRank`, `pCombStephenson` |

## 3. Current Blockers or Open Questions

### Resolved this session

- `make check` failure due to missing `Rfast`, `RcppDist`, and Suggests --- fixed by `make dependencies` (pak).
- `make check` failure due to CMRSS API change --- fixed by translation in
  `pCombStephenson` (commit `cb3604e`).
- LP-solver missing --- installed `highs` 1.12.0-3 and added to Suggests.

### Still open

- **Uncommitted `alpha_adaptive.R` work** (see section 2). Decide whether to
  commit, revise, or revert. Two new untracked tests describe their intent
  in the filenames (`test_irregular_tree_warning.R`,
  `test_weak_fwer_global_null.R`) but are not mentioned in any prior handoff.
- **Push `cb3604e` to `origin/main`?** Currently 1 commit ahead, unpushed.
- **Roxygen-version mismatch.** `make check` no longer warns (clean), but
  `DESCRIPTION` declares `RoxygenNote: 7.3.3` while installed roxygen2 is
  8.0.0. If `make document` is ever re-run, it will bump `RoxygenNote` and
  regenerate the three modified `.Rd` files. Don't run `make document` until
  the uncommitted `alpha_adaptive.R` changes are reviewed --- regen would
  obscure deliberate edits.
- **`pPolyRank` as default `pfn` in `find_blocks()`?** Currently `pIndepDist`.
  User decision, inherited from prior handoff.
- **`pCombStephenson` speed in tree search.** ~0.5s per call; potentially a
  bottleneck for large trees. No fast alternative for quantile-of-effects
  hypotheses (`k < n`) yet. Inherited.
- **`REFACTOR_PLAN.md` cleanup.** Still in repo root (`.Rbuildignore`-ed).
  Inherited.

## 4. Important Context to Preserve

### CMRSS 0.2.6+ k semantics (supersedes any older "footgun" note)

The user-facing `k` argument of `pCombStephenson` is paper-notation
(1..n, indexing the rank of `tau` among all units). The CMRSS internal
argument is now 1..m where `m = sum(Z)`. The wrapper computes
`cmrss_k = k - (n - m)` before calling CMRSS. The hypothesis being tested
is unchanged; only the parameterization moved.

Concrete numerics with `idat` (n = 1000, m = 500):

- Sharp null: user passes `k = 1000` (or omits) → `cmrss_k = 500 = m` → no
  treated unit gets Inf → standard stratified rank-sum.
- Degenerate: user passes any `k <= 500` → wrapper warns, returns `p = 1`,
  never calls CMRSS.
- Quantile-of-effects (`k = 501..999`) → `cmrss_k = 1..499` → LP solves the
  optimal set of `m - cmrss_k` treated units to mask.

### LP solver is now mandatory for non-degenerate CMRSS calls

CMRSS 0.2.6+ calls `solve_optimization()` unconditionally inside
`pval_comb_block`, even when `p = 0`. So *every* non-degenerate call now
needs `highs` or `gurobi`. The 4 CMRSS tests that hit the LP path are
guarded with `skip_if_not_installed("highs")`.

### Makefile structure

The Makefile uses a pattern where targets set `FUNC` (and optionally
`DEVTOOLSARG`) and depend on `.devtools`, which executes
`devtools:::$(FUNC)($(DEVTOOLSARG))`. Adding a new devtools-backed target:

```makefile
.PHONY: newtarget
newtarget: FUNC=devtools_function_name
newtarget: DEVTOOLSARG=arg1 = 'value', arg2 = TRUE
test check check-fast check-cran document build newtarget: .devtools
```

The `dependencies` target is the lone outlier --- it runs pak directly
rather than routing through `.devtools`.

### Quarto vignettes

`DESCRIPTION` declares `VignetteBuilder: quarto`. The Quarto CLI is at
`/usr/local/bin/quarto` (v1.9.37). Sources in `vignettes/`:
`getting-started.qmd`, `advanced-methodologies.qmd`,
`hierarchical-testing-workflow.qmd`. On a machine without Quarto, use
`make check-fast`.

### Test data

`make_test_data.R` creates `idat` with `n = 1000`, `m = 500`, and three
outcomes:

- `Y`: heterogeneous/canceling effects (positive in some blocks, negative
  in others)
- `Yhomog`: strong uniform positive effects across all blocks
- `Ynull`: no treatment effects

### Prior-session context (preserved)

#### Dependency lightening (v0.0.4.1003)

Eight packages moved from Imports to Suggests, guarded with
`requireNamespace()`. `ClusterR` removed entirely. `dataPreparation` replaced
with internal `which_are_constant()` in `R/utils.R`. Imports is now only
`Rcpp`, `Rfast`. Suggests now also includes `CMRSS` and (as of `cb3604e`)
`highs`.

#### Polynomial rank functions

- `pPolyRank()`: `coin::independence_test()` with multivariate polynomial
  rank scores (r=2, 6, 10). Fast (~0.02s). Tests sharp null. Follows
  `pIndepDist` pattern. 25x faster than `pCombStephenson`.
- `pCombStephenson()`: wraps `CMRSS::pval_comb_block()` with the
  `cmrss_k = k - (n - m)` translation. ~0.5s per call.
- Polynomial scores: `(rank(Y) / (n + 1))^(r - 1)`. r = 2 is Wilcoxon-like;
  larger r emphasizes upper-tail effects.

#### Budget weights and switching corollary

`R/alpha_adaptive.R` supports budget-weighted alpha allocation
(`budget_weights` parameter in `compute_adaptive_alphas_tree` and
`alpha_adaptive_tree_pruned`) and a switching corollary (`switching = TRUE`
in `alpha_adaptive_tree_pruned`) that reverts to nominal alpha when the
pruned error load fits within the remaining budget. Tests in
`test_budget_weights.R` and `test_switching_corollary.R`.

#### prune_tree infinite-loop pattern

The `prune_tree()` helper in `tests/testthat/test_switching_corollary.R` had
a bug where descendants were re-discovered every iteration. Fixed by adding
`& !(node_dat$nodenum %in% to_remove)`. Test-only helper.

### Development workflow (CLAUDE_CODING.md)

Checkpoints: (1) write tests first, pause for review; (2) after
implementation, pause before running check; (3) after passing check, pause
for review.

## 5. What's Done vs. What Remains

### Done this session

- Switched `make dependencies` to pak via `pak::pkg_install('.',
  dependencies = TRUE, upgrade = FALSE, ask = FALSE)`.
- Added `make check-fast` and `make check-cran` targets.
- Ran `make dependencies` end-to-end and installed all CRAN deps, the three
  GitHub Remotes (`dsrobertson/onlineFDR`, `bowers-illinois-edu/TreeTestSim`,
  `bowers-illinois-edu/CMRSS`), and `Rfast`/`RcppDist`.
- Diagnosed 6 test failures in `test_stephenson.R` as a CMRSS 0.2.6+ API
  change. Confirmed by reading CMRSS source.
- Patched `pCombStephenson` with `cmrss_k = k - (n - m)` translation,
  validation, and short-circuit return for degenerate `k`.
- Updated `pCombStephenson` docstring.
- Updated `test_stephenson.R` direct-CMRSS comparison and added
  `skip_if_not_installed("highs")` guards.
- Installed `highs` 1.12.0-3 (open-source LP solver) via pak.
- Added `highs` to `DESCRIPTION` Suggests.
- Added `# manytestsr 0.0.4.1005` section to `NEWS.md`.
- Updated `MEMORY.md` `pCombStephenson` description.
- Committed all session changes as `cb3604e`. Pre-existing alpha-adaptive
  work was left uncommitted.
- Verified `make check` passes clean: **0 errors, 0 warnings, 0 notes,
  ~1m 37s, vignettes built, 920 + 10 tests passing, 10 expected skips.**

### Remains

1. **Decide what to do about the uncommitted `alpha_adaptive.R` changes and
   the two untracked test files.** Read the diffs first; the prior handoff
   incorrectly claimed the tree was clean, so the intent of these changes
   is undocumented. They pass `make check`, but that doesn't tell us
   whether they reflect the user's final intent.
2. **Push `cb3604e` to `origin/main`** if and when ready.
3. **Consider running `make document`** once the alpha-adaptive review is
   done, to sync `RoxygenNote:` and any stale `.Rd` files. Don't do this
   before reviewing the uncommitted `.Rd` changes.
4. **Decide on `pPolyRank` as default `pfn`** in `find_blocks()`.
5. **Delete `REFACTOR_PLAN.md`** if the refactor is confirmed complete.
6. **Optional: address `pCombStephenson` speed** as tree-search bottleneck
   (no concrete plan yet; would need a fast alternative for
   quantile-of-effects hypotheses).
