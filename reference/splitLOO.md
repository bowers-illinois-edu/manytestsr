# Splitting function: Leave One Out

A splitting function takes block ids and block ordering vector (or
vectors) and produces a factor that assigns some block ids to one group
or another group.

## Usage

``` r
splitLOO(bid, x)
```

## Arguments

- bid:

  Block id

- x:

  A vector that we can use to order the blocks. A block-level value.
  Like N or harmonic mean weight of the block.

## Value

A factor categorizing blocks into groups.

## Examples

``` r
# Leave-one-out splitting - focuses on largest block vs rest
block_ids <- c("B1", "B2", "B3", "B4", "B5")
block_weights <- c(0.1, 0.3, 0.8, 0.2, 0.4) # B3 is largest

groups <- splitLOO(block_ids, block_weights)

# Show which block is isolated (should be B3 with weight 0.8)
data.frame(block_id = block_ids, weight = block_weights, group = groups)
#>   block_id weight group
#> 1       B1    0.1     1
#> 2       B2    0.3     1
#> 3       B3    0.8     0
#> 4       B4    0.2     1
#> 5       B5    0.4     1
```
