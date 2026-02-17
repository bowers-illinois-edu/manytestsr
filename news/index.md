# Changelog

## manytestsr 0.0.4.1002

### New features

- New exported factory function
  [`alpha_adaptive_tree_pruned()`](../reference/alpha_adaptive_tree_pruned.md)
  creates a branch-pruning adaptive alpha system for use with
  [`find_blocks()`](../reference/find_blocks.md). Unlike
  [`alpha_adaptive_tree()`](../reference/alpha_adaptive_tree.md), which
  pre-computes a fixed schedule, this version can recompute the schedule
  on the surviving subtree after each depth — giving more alpha to
  surviving branches when dead branches are removed. Returns a list with
  three components: `$alphafn` (standard closure), `$update` (recompute
  on pruned tree), and `$reset` (restore full tree).

- [`find_blocks()`](../reference/find_blocks.md) now supports
  list-valued `alphafn` parameters. When `alphafn` is a list (as
  returned by
  [`alpha_adaptive_tree_pruned()`](../reference/alpha_adaptive_tree_pruned.md)),
  `find_blocks` extracts the `$alphafn`, `$update`, and `$reset`
  components, calls `reset` at the start of each run, and calls `update`
  after each depth’s testable decisions. Plain function `alphafn` values
  continue to work unchanged.

- New internal helper
  [`.get_all_descendants()`](../reference/dot-get_all_descendants.md)
  performs BFS traversal on tree-structured `node_dat` to find all
  descendants of given nodes.

## manytestsr 0.0.4.1001

### New features

- New exported function
  [`compute_adaptive_alphas_tree()`](../reference/compute_adaptive_alphas_tree.md)
  computes per-depth adjusted significance levels from an actual
  (possibly irregular) tree structure. Takes `node_dat` with per-node
  sample sizes (as returned by
  [`find_blocks()`](../reference/find_blocks.md)) instead of assuming a
  regular k-ary tree. The algorithm divides alpha at each depth by the
  sum of path powers — the expected number of tests conducted at that
  depth. For regular k-ary trees, this produces identical results to the
  parametric
  [`compute_adaptive_alphas()`](../reference/compute_adaptive_alphas.md).

- New exported factory function
  [`alpha_adaptive_tree()`](../reference/alpha_adaptive_tree.md) creates
  a closure for use with `find_blocks(alphafn = ...)`, using the
  tree-based alpha schedule from
  [`compute_adaptive_alphas_tree()`](../reference/compute_adaptive_alphas_tree.md).
  Drop-in replacement for
  [`alpha_adaptive()`](../reference/alpha_adaptive.md) when the tree has
  irregular branching or unequal sample sizes across nodes.

### Bug fixes

- Fixed 11 test failures in `test_alpha_adaptive.R` that referenced the
  removed `tau` parameter. Two tests (“tau = 0” and “tau = 1”) were
  rewritten as error-load equivalents; the rest had `tau` arguments
  removed.

## manytestsr 0.0.4.1000

### New features

- New exported function
  [`compute_error_load()`](../reference/compute_error_load.md) computes
  the error load at each tree level — the expected number of all-null
  sibling groups that the procedure tests. When the total error load is
  at most 1, the unadjusted procedure controls FWER via natural gating;
  when it exceeds 1, adaptive alpha adjustment is required. Supports
  both a parametric interface (regular k-ary trees with equal splits)
  and a tree interface (irregular trees with per-node sample sizes from
  [`find_blocks()`](../reference/find_blocks.md)).

### Changes

- [`compute_adaptive_alphas()`](../reference/compute_adaptive_alphas.md)
  now checks the error load before computing adjusted alphas. When the
  total error load is at most 1 (natural gating suffices), nominal alpha
  is returned at every level without adjustment. The `tau` parameter has
  been removed; the error load check replaces it.

- [`compute_adaptive_alphas()`](../reference/compute_adaptive_alphas.md)
  gains an `"error_load"` attribute on its return value, so callers can
  inspect the error load diagnostics without a separate call to
  [`compute_error_load()`](../reference/compute_error_load.md).

## manytestsr 0.0.4.0000

### New features

- Added adaptive alpha adjustment for tree-structured hypothesis testing
  ([`compute_adaptive_alphas()`](../reference/compute_adaptive_alphas.md)
  and [`alpha_adaptive()`](../reference/alpha_adaptive.md)). This
  implements Algorithm 1 from the paper’s Appendix D, which adjusts
  significance levels at each tree depth based on estimated power decay.
  When cumulative power is high, alpha is tightened to account for the
  multiplicity of tests that power enables. When cumulative power drops
  below a threshold, natural gating suffices and nominal alpha is used.
  The procedure supports both constant and variable branching factors.

### Interface changes

- The `alphafn` interface used by
  [`find_blocks()`](../reference/find_blocks.md) now passes a `depth`
  parameter (integer vector of tree depths, 1 = root) to alpha
  adjustment functions. This enables alpha strategies that depend on
  tree structure rather than treating p-values as a flat stream.

- [`alpha_investing()`](../reference/alpha_investing.md),
  [`alpha_saffron()`](../reference/alpha_saffron.md), and
  [`alpha_addis()`](../reference/alpha_addis.md) now accept a `depth`
  argument for interface compatibility. They do not use it; their
  behavior is unchanged.

## manytestsr 0.0.3.0000

- Initial tracked version with
  [`find_blocks()`](../reference/find_blocks.md), splitting functions
  (`splitCluster`, `splitEqualApprox`, `splitLOO`,
  `splitSpecifiedFactor`, `splitSpecifiedFactorMulti`), p-value
  functions (`pOneway`, `pWilcox`, `pIndepDist`, `pCombCauchyDist`,
  `pTestTwice`), online FDR alpha adjustment (`alpha_investing`,
  `alpha_saffron`, `alpha_addis`), local p-value adjustment
  (`local_hommel_all_ps`, `local_simes`, `local_bh_all_ps`), and
  reporting/visualization functions.
