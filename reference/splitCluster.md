# Splitting function: K-Means Clustering

A splitting function takes block ids and a block ordering vector (or
vectors) and produces a factor that assigns some block ids to one group
or another group.

## Usage

``` r
splitCluster(bid, x)
```

## Arguments

- bid:

  Block id

- x:

  A vector that we can use to order the blocks

## Value

A factor categorizing blocks into groups.

## Examples

``` r
# Simple example with block weights
block_ids <- c("B1", "B2", "B3", "B4", "B5", "B6")
block_weights <- c(0.1, 0.8, 0.2, 0.9, 0.3, 0.7)

# Split blocks into clusters based on weights
groups <- splitCluster(block_ids, block_weights)
print(groups)
#> [1] 1 0 1 0 1 0
#> Levels: 0 1

# View which blocks are in each group
data.frame(block_id = block_ids, weight = block_weights, group = groups)
#>   block_id weight group
#> 1       B1    0.1     1
#> 2       B2    0.8     0
#> 3       B3    0.2     1
#> 4       B4    0.9     0
#> 5       B5    0.3     1
#> 6       B6    0.7     0
```
