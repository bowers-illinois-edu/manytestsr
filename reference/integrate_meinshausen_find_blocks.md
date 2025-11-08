# Integration function for find_blocks with Meinshausen testing

Integrates Meinshausen hierarchical testing into the find_blocks
workflow

## Usage

``` r
integrate_meinshausen_find_blocks(
  node_dat,
  node_tracker,
  alpha = 0.05,
  method = "simes",
  use_sequential = TRUE
)
```

## Arguments

- node_dat:

  Node-level data from find_blocks

- node_tracker:

  Node tracking structure

- alpha:

  Error rate

- method:

  P-value combination method

- use_sequential:

  Use sequential rejection principle

## Value

Enhanced node_dat with Meinshausen results
