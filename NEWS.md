# manytestsr 0.0.4.0000

## New features

* Added adaptive alpha adjustment for tree-structured hypothesis testing
  (`compute_adaptive_alphas()` and `alpha_adaptive()`). This implements
  Algorithm 1 from the paper's Appendix D, which adjusts significance
  levels at each tree depth based on estimated power decay. When
  cumulative power is high, alpha is tightened to account for the
  multiplicity of tests that power enables. When cumulative power drops
  below a threshold, natural gating suffices and nominal alpha is used.
  The procedure supports both constant and variable branching factors.

## Interface changes

* The `alphafn` interface used by `find_blocks()` now passes a `depth`
  parameter (integer vector of tree depths, 1 = root) to alpha
  adjustment functions. This enables alpha strategies that depend on
  tree structure rather than treating p-values as a flat stream.

* `alpha_investing()`, `alpha_saffron()`, and `alpha_addis()` now accept
  a `depth` argument for interface compatibility. They do not use it;
  their behavior is unchanged.

# manytestsr 0.0.3.0000

* Initial tracked version with `find_blocks()`, splitting functions
  (`splitCluster`, `splitEqualApprox`, `splitLOO`, `splitSpecifiedFactor`,
  `splitSpecifiedFactorMulti`), p-value functions (`pOneway`, `pWilcox`,
  `pIndepDist`, `pCombCauchyDist`, `pTestTwice`), online FDR alpha
  adjustment (`alpha_investing`, `alpha_saffron`, `alpha_addis`), local
  p-value adjustment (`local_hommel_all_ps`, `local_simes`,
  `local_bh_all_ps`), and reporting/visualization functions.
