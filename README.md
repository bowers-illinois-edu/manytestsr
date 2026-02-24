# manytestsr

Testing to detect and localize heterogeneous effects in block-randomized experiments.

This package implements the nesting (top-down) procedures detailed in
[Bowers and Chen (2020)](https://doi.org/10.1093/pan/mpaa031). It recursively
splits experimental blocks, tests for treatment effects at each node, and
controls familywise error using adaptive alpha adjustment — so you can ask
not just *whether* a treatment worked, but *where* it worked.

## Warning: development-stage package

**This package is under active development. The API may change without
notice, and not all procedures have proven statistical properties yet.** If
you use it in applied work, pin a specific commit or version and check back
for updates.

Current version: **0.0.4.1002**

## Installation

```r
# install.packages("remotes")
remotes::install_github("bowers-illinois-edu/manytestsr")
```

The package includes C++ code (via Rcpp/RcppArmadillo) and requires a working
C++ compiler. On macOS, install Xcode Command Line Tools; on Windows, install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/).

## Public TODO list

Items marked with a check are done. Items without a check are open — they
represent known limitations or planned work.

### Alpha adjustment and error control

- [x] Adaptive alpha adjustment for tree-structured testing
      (`compute_adaptive_alphas`, `alpha_adaptive`)
- [x] Tree-based alpha schedules for irregular trees
      (`compute_adaptive_alphas_tree`, `alpha_adaptive_tree`)
- [x] Branch-pruning adaptive alpha that reallocates alpha from dead branches
      (`alpha_adaptive_tree_pruned`)
- [x] Error-load diagnostics (`compute_error_load`)
- [ ] Remove or deprecate the sequential (stream-based) alpha procedures
      (`alpha_investing`, `alpha_saffron`, `alpha_addis`) from `onlineFDR`.
      These treat p-values as a flat stream and do not have provable FWER
      control properties for tree-structured testing.

### Local p-value adjustment

- [ ] Modify the "local Hommel" style adjustments (`local_hommel_all_ps`) to
      operate level-by-level rather than node-by-node.

### Testing

- [ ] Incorporate a test of the sharp null of no effect for any unit, from the
      [CMRSS](https://github.com/bowers-illinois-edu/CMRSS) package.

### Documentation and usability

- [ ] Complete the vignettes (`getting-started`, `hierarchical-testing-workflow`,
      `advanced-methodologies`).
- [ ] Add a worked example to this README.

## Key functions

| Function | Purpose |
|----------|---------|
| `find_blocks()` | Core recursive splitting and testing procedure |
| `compute_adaptive_alphas()` | Depth-adjusted significance levels for regular trees |
| `alpha_adaptive_tree_pruned()` | Branch-pruning alpha for irregular trees |
| `compute_error_load()` | Diagnose whether natural gating controls FWER |
| `splitCluster`, `splitEqualApprox`, `splitLOO`, `splitSpecifiedFactor` | Splitting strategies |
| `pOneway`, `pWilcox`, `pIndepDist`, `pCombCauchyDist` | P-value functions for different test statistics |
| `local_hommel_all_ps`, `local_simes`, `local_bh_all_ps` | Local p-value adjustment |
| `make_results_tree`, `make_results_ggraph` | Visualization of results |

## Development

```
make document   # generate roxygen2 documentation
make test       # run test suite
make check      # R CMD check
make build      # build the package
```

The project uses [renv](https://rstudio.github.io/renv/) for dependency
management. Run `make dependencies` to install required packages.

## Implementation notes

Distance and transformation calculations use C++ (Rcpp/RcppArmadillo) with
three code paths selected by dataset size:

- `fast_dists_by_unit_arma2_par` — OpenMP parallel processing (`parallel = "yes"`)
- `fast_dists_and_trans` — direct matrix computation for small N
- `fast_dists_and_trans_by_unit_arma` — unit-by-unit computation for N > 20,
  avoiding large matrices in memory
