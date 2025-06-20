#' @useDynLib manytestsr, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.formula dist kmeans median na.omit p.adjust rnorm runif sd terms var
NULL

# a fix for R check warnings from https://github.com/STAT545-UBC/Discussion/issues/451
# see for example https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

utils::globalVariables(c(
  "J", "Y", "a", "allnull", "alpha1", "anynotnull", "bF", "batch", "biggrp", "blocksbygroup", "both_ns",
  "depth", "fin_grp", "fin_nodenum", "fin_parent", "fin_parent_p", "g1", "group_hit",
  "group_hit2", "grpsize", "hit", "hit_grp", "i.nodesize", "label", "leaves", "max_alpha",
  "max_p", "maxate", "maxdepth", "meanate", "medianate", "minate", "name", "newZ", "newZF",
  "nodenum", "nodenum_current", "nodenum_prev", "nodes", "nodesize", "out_degree", "p", "p1",
  "parent_alpha", "parent_name", "parent_of_ns", "paths", "pfinalb", "shortbf",
  "single_hit", "testable", "trueate", "trueblocks", "truetaui", "vecdist", "y1new", "y1sim"
))

## `utils` is loaded first to ensure proper ordering in the search
## stack, so that devtool's versions of `help` and `?` mask `utils`.
