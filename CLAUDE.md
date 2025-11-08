# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

This is an R package called `manytestsr` that implements nesting or
top-down procedures for testing to detect heterogeneous effects in
block-randomized experiments. The package combines statistical testing
with adaptive splitting procedures to localize causal effects.

## Development Commands

### Core Development Workflow

- `make document` - Generate documentation using roxygen2
- `make test` - Run test suite (excludes profiling tests)
- `make check` - Run R CMD check for package validation
- `make build` - Build the package
- `make dependencies` - Install package dependencies

### Interactive Development

- `make interactive` - Start R session with package loaded via `load.R`
- `R -q --no-save` with `R_PROFILE=load.R` for development session

### Testing

- Tests are in `tests/testthat/` directory
- Main test runner: `tests/testthat.R` uses
  `test_check("manytestsr", filter = "profil", invert = TRUE)`
- Set `interactive <- TRUE` in test files for debugging (should be
  `FALSE` for production)
- Tests use `setDTthreads(1)` to ensure reproducible results

### Build System

- Uses `renv` for dependency management (see `renv.lock`)
- C++ compilation via Rcpp/RcppArmadillo with OpenMP support
- Makefile provides convenient development targets

## Code Architecture

### Core Components

1.  **Statistical Testing Functions** (`R/pval_fns.R`)
    - P-value computation functions for different test statistics
    - Permutation and asymptotic testing approaches
2.  **Distance and Transformation Functions** (`R/dists.R`,
    `src/dists_and_trans.cpp`)
    - Fast C++ implementations for computing outcome distances
    - Three main C++ functions based on data size:
      - `fast_dists_by_unit_arma2_par` for parallel processing (OpenMP)
      - `fast_dists_and_trans` for small n
      - `fast_dists_and_trans_by_unit_arma` for N\>20
    - Includes energy distance calculations and rank-based
      transformations
3.  **Block Finding and Splitting** (`R/find_blocks.R`,
    `R/splitting_fns.R`)
    - Core [`find_blocks()`](reference/find_blocks.md) function
      implements the adaptive testing procedure
    - Various splitting strategies (equal, LOO, specified,
      cluster-based)
    - Tree-based approach with block identification and testing
4.  **Alpha Adaptation** (`R/alpha_fns.R`)
    - Implements online FDR control procedures
    - Wrappers for `onlineFDR` package functions (ADDIS, Saffron,
      alpha-investing)
5.  **P-value Adjustment** (`R/p_adjustment.R`)
    - Local adjustment procedures (Simes, Hommel, min-p, etc.)
    - Integration with global alpha spending procedures

### Key Dependencies

- **Core**: `data.table`, `Rcpp`, `RcppArmadillo` for performance
- **Statistical**: `coin` for exact tests, `hommel` for p-value
  adjustment
- **Visualization**: `ggplot2`, `ggraph`, `tidygraph` for result
  plotting
- **Clustering**: `ClusterR`, `Ckmeans.1d.dp` for adaptive splitting

### Data Flow

1.  Input data at unit (`idat`) and block (`bdat`) levels
2.  Recursive splitting using various strategies guided by `splitby`
    variable
3.  Statistical testing at each node using permutation or asymptotic
    methods
4.  Alpha adaptation using online FDR procedures
5.  Local p-value adjustment within nodes
6.  Tree construction with results stored for reporting

### Performance Considerations

- C++ implementations for computationally intensive distance
  calculations
- OpenMP parallelization available for large datasets
- Memory-efficient approaches for large N (avoids holding large
  matrices)
- Uses `data.table` for efficient data manipulation
- Thread control via `setDTthreads(1)` for reproducible results

### Testing Strategy

- Comprehensive test suite covering all major functions
- Tests use simulated data from `TreeTestSim` package
- Performance profiling code available in `tests/` directory
- Tests exclude profiling by default
  (`filter = "profil", invert = TRUE`)
