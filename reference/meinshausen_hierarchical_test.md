# Meinshausen's Hierarchical Testing with Sequential Rejection Principle

Implementation of Meinshausen's (2008) hierarchical testing of variable
importance enhanced with the sequential rejection principle of Goeman
and Solari (2010). This approach provides a unified framework for
hierarchical testing that maintains FWER control while leveraging the
hierarchical structure for increased power.

## Usage

``` r
meinshausen_hierarchical_test(
  node_dat,
  node_tracker,
  alpha = 0.05,
  method = "simes",
  use_sequential = TRUE
)
```

## Arguments

- node_dat:

  Data.table containing node-level test results with columns: nodenum,
  p, depth, nodesize, and optionally parent, testable

- node_tracker:

  Node tracking object containing the hierarchical structure

- alpha:

  Global Type I error rate (default: 0.05)

- method:

  Method for combining p-values ("simes", "fisher", "bonferroni")

- use_sequential:

  Logical, whether to use sequential rejection principle

## Value

Updated node_dat with Meinshausen hierarchical testing results

## References

Meinshausen, N. (2008). Hierarchical testing of variable importance.
Biometrika 95, 265-278.

Goeman, J. J., & Solari, A. (2010). The sequential rejection principle
of familywise error control. Annals of Statistics 38, 3782-3810.
