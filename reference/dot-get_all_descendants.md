# Find All Descendants of Given Nodes in a Tree

BFS traversal using the parent column to collect every node descended
from the given starting nodes.

## Usage

``` r
.get_all_descendants(node_dat, nodenums)
```

## Arguments

- node_dat:

  data.frame with at least `nodenum` and `parent` columns.

- nodenums:

  Integer vector of node IDs whose descendants to find.

## Value

Integer vector of descendant node IDs (not including the starting nodes
themselves). Empty integer vector if no descendants exist.
