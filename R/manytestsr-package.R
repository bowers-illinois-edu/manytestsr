#' @useDynLib manytestsr, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.formula dist kmeans median na.omit p.adjust rnorm runif sd terms var as.dist dt hclust pbinom pchisq pnorm pt setNames
#' @importFrom utils combn
NULL

# a fix for R check warnings from https://github.com/STAT545-UBC/Discussion/issues/451
# see for example https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

utils::globalVariables(c(
  "J", "Y", "a", "allnull", "alpha1", "anynotnull", "bF", "batch", "group_id", "block_id", "blocks", "blocksbygroup", "both_ns",
  "closed_testing_reject", "consistent", "dat", "decision.times", "depth", "desc_min_p", "dt",
  "fin_grp", "fin_nodenum", "fin_parent", "fin_parent_p", "from", "g1", "gamma", "group_hit",
  "group_hit2", "grpsize", "hit", "hit_grp", "i.blocks", "i.node_label_val", "i.nonnull", "i.nodesize", "i.num_leaves", "id",
  "iu_p_intersection", "iu_p_union", "iu_reject_intersection", "iu_reject_union",
  "label", "leaves", "max_alpha", "max_p", "maxate", "maxdepth", "meanate", "medianate", "meinshausen_adjusted_p",
  "meinshausen_level_tested", "meinshausen_reject", "minate", "name", "nb", "newZ", "newZF",
  "node_id", "node_number", "node_type", "nodenum", "nodenum_current", "nodenum_prev", "nodes", "nodesize", "nonnull",
  "object", "original_name", "out_degree", "p", "p1", "p_intersection", "p_union", "p_value_lower", "p_value_upper",
  "parent", "parent_alpha", "parent_id", "parent_name", "parent_of_ns", "paths", "pfinalb", "rankY",
  "reject_intersection", "reject_union", "shortbf", "single_hit", "testable", "thetreatD", "to",
  "trueate", "trueblocks", "truetaui", "variables", "vecdist", "y1new", "y1sim"
))

## `utils` is loaded first to ensure proper ordering in the search
## stack, so that devtool's versions of `help` and `?` mask `utils`.
