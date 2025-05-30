% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reporting.R
\name{make_results_tree}
\alias{make_results_tree}
\title{Make a node-level dataset from a block-level dataset}
\usage{
make_results_tree(
  orig_res,
  block_id = NULL,
  node_label = NULL,
  return_what = "all"
)
}
\arguments{
\item{orig_res}{data.table from find_blocks(); must include elements such as
\itemize{
\item biggrp    (dot-sep lineage, may be truncated),
\item p1,p2,… and a pfinal*,
\item alpha1, alpha2, …
}}

\item{block_id}{optional name of your block ID column (e.g. "bF"), or NULL}

\item{node_label}{optional name of a descriptive label column}

\item{return_what}{a character vector containing "all", "graph" (a tbl_graph
object with nodes and edges), "nodes" (a data.table with node level
information), "test_summary" (a data.table object with one row indicating
false and true discoveries, etc.)}
}
\value{
a list that can contain nodes, a tbl_graph object, and/or a test_summary
}
\description{
Given the results of the splitting and testing algorithm, make a node level
data set for use in reporting results and as input to ggraph for
visualization in terms of a tree graph.
}
