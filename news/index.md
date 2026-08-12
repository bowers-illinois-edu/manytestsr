# Changelog

## manytestsr 0.0.4.1008

Corrections responding to an anonymous AOAS referee who worked through
this package alongside the manuscript and identified genuine errors
(2026-08; see FIX_PLAN.md for the mapping). The pre-correction state is
tagged `pre-referee-fixes` so previously published numbers remain
reproducible.

### Breaking changes

- [`compute_error_load()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_error_load.md)
  and the adaptive-alpha schedule builders now implement the paper’s
  Definition 2: per-depth loads sum path power (product of
  PROPER-ancestor rejection probabilities) over the nodes at each depth,
  from depth 2. The previous formula multiplied each node’s own theta
  back in and included a depth-1 term, understating the load and
  returning `needs_adjustment = FALSE` in designs that require
  adjustment. The `needs_adjustment` gate and the alpha-schedule
  denominators now use the same quantity.
- Power calculations behind the error load are now two-tailed:
  `theta = pnorm(delta*sqrt(n) - z) + pnorm(-delta*sqrt(n) - z)`, so
  theta at `delta = 0` equals the size alpha rather than alpha/2.
- The error-load/power path refuses `nodesize` values below 1 with an
  error.
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)’s
  default `blocksize = "hwt"` stores harmonic weights whose node sums
  lie on the unit interval; feeding them to a power formula produced
  meaningless near-floor thetas. Pass a headcount column
  (e.g. `blocksize = "nb"`) for any error-load or adaptive-alpha use;
  `hwt` remains valid for splitting and weighting.
- [`alpha_adaptive_tree_pruned()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree_pruned.md)
  now warns at creation that the strong-FWER guarantee for the
  pruned-load schedule and its switching rule DOES NOT HOLD: the
  underlying theorem was falsified by exact counterexample (FWER 0.063
  at alpha 0.05 with all hypotheses satisfied). Levels are computed as
  before for reproducibility; for a schedule with a proof use
  [`alpha_adaptive_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree.md)
  with static budget weights.

### Bug fixes

- [`report_detections()`](https://bowers-illinois-edu.github.io/manytestsr/reference/report_detections.md)
  computed each block’s parent p-value at the GLOBAL maximum depth, so
  branches whose testing stopped earlier got the wrong parent or NA
  (surfacing as `hit = NA`, silently dropped by `sum(hit, na.rm = TRUE)`
  callers), and group coverage was missed on any branch shorter than the
  deepest one. All quantities now use each block’s own final tested
  depth, and `hit` is never NA.
- [`report_detections()`](https://bowers-illinois-edu.github.io/manytestsr/reference/report_detections.md)
  gains `hit_type` (“single” = the block’s own test rejected; “group” =
  family parent rejected with no child rejected; “none” otherwise – a
  sibling’s individual rejection explains the parent, so coverage does
  not spread to failed siblings) and `group_p` (the rejecting parent’s
  p-value for covered blocks). `pfinalb` is now the running maximum of
  p-values along the block’s tested path rather than the deepest p
  reached.
- Removed a stray positional argument in five
  [`coin::approximate()`](https://rdrr.io/pkg/coin/man/NullDistribution.html)
  calls that bound to coin’s `cl` formal; harmless under
  `parallel = "no"` but broken under `parallel = "snow"`.

## manytestsr 0.0.4.1007

### Bug fixes

- [`make_results_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/make_results_tree.md)
  now propagates block-level truth (`nonnull`) up the tree instead of
  labeling leaf nodes only. An internal node is non-null iff a
  descendant leaf is non-null, and a known null iff all of its
  descendant leaves are known nulls. Without this, a false rejection of
  a true internal null hypothesis (rejecting a whole null
  group/college/region) was dropped from the node-level FWER tally
  (`node_any_false_rejection`, `node_false_rejection_prop`,
  `node_num_false_rejections`, `node_false_discovery_prop`), and correct
  rejections of non-null ancestors were missing from
  `node_true_discoveries` and `node_power`. The procedure’s rejection
  decisions are unchanged; only the truth labeling used to score
  node-level metrics is corrected. New tests in
  `tests/testthat/test_node_truth_propagation.R`.

## manytestsr 0.0.4.1005

### Bug fixes

- [`pCombStephenson()`](https://bowers-illinois-edu.github.io/manytestsr/reference/pCombStephenson.md)
  is now compatible with `CMRSS` 0.2.6+, which constrains its internal
  `k` argument to `1..m` (where `m = sum(Z)` is the number of treated
  units) rather than `1..n`. The wrapper’s user-facing `k` argument
  remains paper-notation in `1..n` (the rank index of tau among all
  units), and the wrapper now translates internally as
  `cmrss_k = k - (n - m)` before calling
  [`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html).
  The default `k = n` still tests Fisher’s sharp null. The degenerate
  path (`k <= n - m`) now warns and returns `p = 1` without calling
  CMRSS, avoiding the LP solver on a branch whose test statistic is
  constant across permutations anyway.

- [`pCombStephenson()`](https://bowers-illinois-edu.github.io/manytestsr/reference/pCombStephenson.md)
  now validates `k <= n` and stops with a clear error if `k` exceeds the
  total number of units.

### Dependency changes

- Added `highs` to Suggests. `CMRSS` 0.2.6+ requires an LP solver
  (`highs` or `gurobi`) for the non-degenerate path of
  `pval_comb_block()`. `highs` is open-source and the recommended
  solver. Install with `install.packages('highs')` or via the
  GitHub-only CMRSS GitHub remote.

## manytestsr 0.0.4.1003

### New features

- New exported function
  [`pPolyRank()`](https://bowers-illinois-edu.github.io/manytestsr/reference/pPolyRank.md)
  tests Fisher’s sharp null of no effects using multiple polynomial rank
  score functions simultaneously via
  [`coin::independence_test()`](https://rdrr.io/pkg/coin/man/IndependenceTest.html).
  Computes within-block polynomial scores at multiple r values (default
  r = 2, 6, 10) and passes them as a multivariate response, providing
  adaptive sensitivity to treatment effects without pre-committing to a
  single rank scoring.

- New exported function
  [`pCombStephenson()`](https://bowers-illinois-edu.github.io/manytestsr/reference/pCombStephenson.md)
  provides a formula-based wrapper around
  [`CMRSS::pval_comb_block()`](https://bowers-illinois-edu.github.io/CMRSS/reference/pval_comb_block.html),
  the combined Stephenson rank test of Kim, Li, and Bowers. Tests
  quantile-of-effects hypotheses (whether the k-th largest individual
  effect exceeds a threshold); defaults to k = n, c = 0 for the sharp
  null. CMRSS is a Suggests dependency installed from GitHub
  (`bowers-illinois-edu/CMRSS`).

### Dependency changes

- Moved 8 packages from Imports to Suggests: `stringi`, `tidygraph`,
  `ggraph`, `digest`, `ggplot2`, `Ckmeans.1d.dp`, `onlineFDR`, `hommel`.
  Core test statistic functions (`pIndepDist`, `pTestTwice`,
  `pCombCauchyDist`, `pOneway`, `pWilcox`) now install without pulling
  in heavy tree-testing and visualization libraries. Functions that need
  the moved packages check with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and give a
  clear error message if the package is missing.

- Removed `ClusterR` from Imports (unused; the code that called it was
  already commented out).

- Replaced `dataPreparation` dependency with an internal
  `which_are_constant()` helper (a one-liner that checks for columns
  with fewer than 2 unique values). `dataPreparation` remains in
  Suggests for cross-validation testing only.

- Added `CMRSS` to Suggests (GitHub-only: `bowers-illinois-edu/CMRSS`).

## manytestsr 0.0.4.1002

### New features

- New exported factory function
  [`alpha_adaptive_tree_pruned()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree_pruned.md)
  creates a branch-pruning adaptive alpha system for use with
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md).
  Unlike
  [`alpha_adaptive_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree.md),
  which pre-computes a fixed schedule, this version can recompute the
  schedule on the surviving subtree after each depth — giving more alpha
  to surviving branches when dead branches are removed. Returns a list
  with three components: `$alphafn` (standard closure), `$update`
  (recompute on pruned tree), and `$reset` (restore full tree).

- [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)
  now supports list-valued `alphafn` parameters. When `alphafn` is a
  list (as returned by
  [`alpha_adaptive_tree_pruned()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree_pruned.md)),
  `find_blocks` extracts the `$alphafn`, `$update`, and `$reset`
  components, calls `reset` at the start of each run, and calls `update`
  after each depth’s testable decisions. Plain function `alphafn` values
  continue to work unchanged.

- New internal helper
  [`.get_all_descendants()`](https://bowers-illinois-edu.github.io/manytestsr/reference/dot-get_all_descendants.md)
  performs BFS traversal on tree-structured `node_dat` to find all
  descendants of given nodes.

## manytestsr 0.0.4.1001

### New features

- New exported function
  [`compute_adaptive_alphas_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas_tree.md)
  computes per-depth adjusted significance levels from an actual
  (possibly irregular) tree structure. Takes `node_dat` with per-node
  sample sizes (as returned by
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md))
  instead of assuming a regular k-ary tree. The algorithm divides alpha
  at each depth by the sum of path powers — the expected number of tests
  conducted at that depth. For regular k-ary trees, this produces
  identical results to the parametric
  [`compute_adaptive_alphas()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas.md).

- New exported factory function
  [`alpha_adaptive_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree.md)
  creates a closure for use with `find_blocks(alphafn = ...)`, using the
  tree-based alpha schedule from
  [`compute_adaptive_alphas_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas_tree.md).
  Drop-in replacement for
  [`alpha_adaptive()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive.md)
  when the tree has irregular branching or unequal sample sizes across
  nodes.

### Bug fixes

- Fixed 11 test failures in `test_alpha_adaptive.R` that referenced the
  removed `tau` parameter. Two tests (“tau = 0” and “tau = 1”) were
  rewritten as error-load equivalents; the rest had `tau` arguments
  removed.

## manytestsr 0.0.4.1000

### New features

- New exported function
  [`compute_error_load()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_error_load.md)
  computes the error load at each tree level — the expected number of
  all-null sibling groups that the procedure tests. When the total error
  load is at most 1, the unadjusted procedure controls FWER via natural
  gating; when it exceeds 1, adaptive alpha adjustment is required.
  Supports both a parametric interface (regular k-ary trees with equal
  splits) and a tree interface (irregular trees with per-node sample
  sizes from
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)).

### Changes

- [`compute_adaptive_alphas()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas.md)
  now checks the error load before computing adjusted alphas. When the
  total error load is at most 1 (natural gating suffices), nominal alpha
  is returned at every level without adjustment. The `tau` parameter has
  been removed; the error load check replaces it.

- [`compute_adaptive_alphas()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas.md)
  gains an `"error_load"` attribute on its return value, so callers can
  inspect the error load diagnostics without a separate call to
  [`compute_error_load()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_error_load.md).

## manytestsr 0.0.4.0000

### New features

- Added adaptive alpha adjustment for tree-structured hypothesis testing
  ([`compute_adaptive_alphas()`](https://bowers-illinois-edu.github.io/manytestsr/reference/compute_adaptive_alphas.md)
  and
  [`alpha_adaptive()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive.md)).
  This implements Algorithm 1 from the paper’s Appendix D, which adjusts
  significance levels at each tree depth based on estimated power decay.
  When cumulative power is high, alpha is tightened to account for the
  multiplicity of tests that power enables. When cumulative power drops
  below a threshold, natural gating suffices and nominal alpha is used.
  The procedure supports both constant and variable branching factors.

### Interface changes

- The `alphafn` interface used by
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)
  now passes a `depth` parameter (integer vector of tree depths, 1 =
  root) to alpha adjustment functions. This enables alpha strategies
  that depend on tree structure rather than treating p-values as a flat
  stream.

- [`alpha_investing()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_investing.md),
  [`alpha_saffron()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_saffron.md),
  and
  [`alpha_addis()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_addis.md)
  now accept a `depth` argument for interface compatibility. They do not
  use it; their behavior is unchanged.

## manytestsr 0.0.3.0000

- Initial tracked version with
  [`find_blocks()`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md),
  splitting functions (`splitCluster`, `splitEqualApprox`, `splitLOO`,
  `splitSpecifiedFactor`, `splitSpecifiedFactorMulti`), p-value
  functions (`pOneway`, `pWilcox`, `pIndepDist`, `pCombCauchyDist`,
  `pTestTwice`), online FDR alpha adjustment (`alpha_investing`,
  `alpha_saffron`, `alpha_addis`), local p-value adjustment
  (`local_hommel_all_ps`, `local_simes`, `local_bh_all_ps`), and
  reporting/visualization functions.
