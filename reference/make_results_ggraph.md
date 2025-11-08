# Make a plot of the nodes

Given the results of the splitting and testing algorithm in the form of
a graph from [make_results_tree](make_results_tree.md), make a node
level data set for use in reporting results in terms of a binary tree
graph. This does not print or plot the graph. You'll need to do that
with the resulting object.

## Usage

``` r
make_results_ggraph(res_graph, remove_na_p = TRUE)
```

## Arguments

- res_graph:

  A tidygraph object produced from make_results_tree

- remove_na_p:

  A logical indicating whether the graph should include nodes/leaves
  that were not tested. Default (TRUE) is to remove them. When
  remove_na_p is FALSE, the graph might look strange since some blocks
  will not have a known position in the graph (the graph is fixed, but
  not specified by the find_blocks function when a node or block is not
  visited for testing.)

## Value

A ggraph object
