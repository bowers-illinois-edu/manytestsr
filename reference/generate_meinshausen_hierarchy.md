# Generate Meinshausen hierarchical clustering for variables

Creates a hierarchical clustering structure suitable for Meinshausen
testing when applied to variable selection problems. This is typically
used when the goal is to test groups of correlated variables
hierarchically.

## Usage

``` r
generate_meinshausen_hierarchy(
  correlation_matrix,
  method = "complete",
  min_cluster_size = 2
)
```

## Arguments

- correlation_matrix:

  Correlation matrix of variables

- method:

  Clustering method ("complete", "average", "single")

- min_cluster_size:

  Minimum size for clusters to be tested

## Value

List with clustering structure compatible with node_tracker format
