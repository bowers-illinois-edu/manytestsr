# Apply consonance correction to test results

Apply consonance correction to test results

## Usage

``` r
apply_consonance_correction(
  node_dat,
  tracker,
  rejection_column = "testable",
  correction_method = "propagate_up"
)
```

## Arguments

- node_dat:

  Node data table with test results

- tracker:

  Node tracker object

- rejection_column:

  Column with rejection decisions

- correction_method:

  Method for correction ("propagate_up", "propagate_down",
  "closed_testing")

## Value

Corrected node data table
