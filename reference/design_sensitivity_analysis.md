# Design Sensitivity Analysis for Hierarchical Testing

Implementation of Rosenbaum's design sensitivity analysis to assess
robustness of treatment effect findings to unobserved confounding. This
extends the standard single-study approach to hierarchical testing.

## Usage

``` r
design_sensitivity_analysis(
  idat,
  bdat,
  blockid,
  formula,
  gamma_range = seq(1, 3, by = 0.1),
  test_statistic = "wilcoxon",
  alpha = 0.05
)
```

## Arguments

- idat:

  Individual-level data

- bdat:

  Block-level data

- blockid:

  Block identifier

- formula:

  Test formula

- gamma_range:

  Vector of sensitivity parameters to test

- test_statistic:

  Type of test statistic ("wilcoxon", "sign", "hodges_lehmann")

- alpha:

  Type I error rate

## Value

Design sensitivity analysis results

## References

Rosenbaum, P. R. (2017). Observation and experiment: an introduction to
causal inference. Harvard University Press.
