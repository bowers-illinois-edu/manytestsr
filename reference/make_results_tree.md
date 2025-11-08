# Make a node-level dataset from a block-level dataset

Given the results of the splitting and testing algorithm, make a node
level data set for use in reporting results and as input to ggraph for
visualization in terms of a tree graph.

## Usage

``` r
make_results_tree(
  orig_res,
  block_id,
  node_label = NULL,
  return_what = "all",
  truevar_name = NULL
)
```

## Arguments

- orig_res:

  data.table from find_blocks(); must include elements such as

  - biggrp (dot-sep lineage, may be truncated),

  - p1,p2,… and a pfinal\*,

  - alpha1, alpha2, …

- block_id:

  name of your block ID column (e.g. "bF")

- node_label:

  optional name of a descriptive label column

- return_what:

  a character vector containing "all", "graph" (a tbl_graph object with
  nodes and edges), "nodes" (a data.table with node level information),
  "test_summary" (a data.table object with one row indicating false and
  true discoveries, etc.)

- truevar_name:

  the optional name of a column recording the true treatment effect
  (used here to find blocks where the true effect is 0 or not). In some
  simulations we have a column called nonnull which is TRUE if that
  block or node has a non-zero effect and FALSE if the block or node has
  a truly zero effect. So, truevar_name can be "nonnull"

## Value

a list that can contain nodes, a tbl_graph object, and/or a test_summary
