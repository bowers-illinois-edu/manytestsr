# Functions for organizing and graphing results

#' Return detected blocks plus info
#'
#' Given the results of the splitting and testing algorithm, report on the
#' blocks where the null of no effects could be rejected at level alpha.
#' Currently calculates rejections using an FWER style criteria (p of a node =
#' max of all previous nodes) if the final alphas are all the same as the
#' scalar alpha OR if fwer=TRUE.

#' @param orig_res results data.table output from the \code{\link{find_blocks}} function.
#' @param fwer (default is TRUE) means that a block is detected (or not) using the maximum p-value associated with the
#' block (or the groups containing that block). fwer=FALSE to detect blocks (or groups of blocks) using FDR control.
#' @param alpha Is the false positive rate used for detecting an effect if it is constant (i.e. not an FDR-style approach).
#' @param only_hits (default FALSE) returns only the detected blocks instead of all of them
#' @param blockid Name of block variable (the blocking variable is a factor)

#' @return A data.table adding a column \code{hit} to the \code{res} data.table
#' indicating a "hit" or detection for that block (or group of blocks)

#' @importFrom stringi stri_count_fixed stri_split_fixed
#' @import data.table
#' @export
report_detections <- function(orig_res, fwer = TRUE, alpha = .05, only_hits = FALSE, blockid = "blockF") {
  res <- copy(orig_res)

  # This next summarizes hits for the output from adjust_block_tests (the bottom-up or test all blocks method)
  if (length(grep("biggrp", names(res))) == 0) {
    # This is for the bottom-up/test every block method
    # max_p are the adjusted p-values so we can use alpha=.05 for error rate control
    res[, hit := max_p <= alpha]
    res[, hit_grp := nodenum_current]
    if (all(!res$hit)) {
      res[, hit_grp := NA]
    }
  } else {
    # For the splitting based methods (output  from find_blocks)
    res[, fin_nodenum := nodenum_current]
    res[, fin_parent := nodenum_prev]
    # Maximum tree depth for a node encoded in the biggrp string: basically number of dots+1 or number of node numbers
    # (where node numbers are separated by a dot)
    res[, maxdepth := stri_count_fixed(biggrp, ".") + 1]
    if (all(res$maxdepth == 1)) {
      # When the algorithm stops at the first level there are no parents.
      # There is just the root node. This is an attempt at a work around
      res[, fin_parent_p := p1]
      res[, parent_alpha := alpha1]
    } else {
      res[, fin_parent_p := get(paste0("p", maxdepth - 1)), by = seq_len(nrow(res))]
      res[, parent_alpha := get(paste0("alpha", maxdepth - 1)), by = seq_len(nrow(res))]
    }

    # If the block is a terminal node and p<=alpha then this is a hit or
    # detected effect If the block is a terminal node and p>alpha for both it
    # and all other terminal nodes then we have a hit or detected effect for
    # *all* terminal nodes but cannot distinguish between them. I think max_p
    # should be p at maxdepth for the FDR algorithm and maximum p overall (or
    # pfinalb) for the FWER algo.

    if (fwer) {
      res[, max_p := pfinalb]
      res[, max_alpha := alpha]
    } else {
      res[, max_p := get(paste0("p", maxdepth)), by = seq_len(nrow(res))]
      res[, max_alpha := get(paste0("alpha", maxdepth)), by = seq_len(nrow(res))]
    }

    # A detection on a single block is scored if the final p <= alpha for a node containing a single block (i.e. a leaf)
    res[, single_hit := (max_p <= max_alpha & blocksbygroup == 1)]

    # A detection is also scored if all descendents have p > alpha but the parent
    # has p <= alpha: this is a grouped detection with multiple blocks.

    ## So, record the minimum p for all descendents of the parents where the splitting has ended
    if (all(res$biggrp == "1")) {
      ## when there is a single group then fin_parent==0 and so doesn't show up in biggrp
      res[, desc_min_p := max_p]
    } else {
      res[, desc_min_p := {
        which_desc <- grep(unique(fin_parent), res$biggrp)
        min(res$max_p[which_desc])
      }, by = fin_parent]
    }

    res[, group_hit := fifelse(!single_hit & (all(max_p > max_alpha) &
      (fin_parent_p <= parent_alpha) & desc_min_p > max_alpha), TRUE, FALSE),
    by = fin_parent
    ]

    res[, hit := single_hit | group_hit]
    if (any(res$hit)) {
      res[(hit), fin_grp := fifelse(single_hit, fin_nodenum, fin_parent)]
      res[(hit), hit_grp := fin_grp]
      res[!(hit), hit_grp := fin_nodenum]
      # Make sure that no  hit_grp includes *both* detections and misses/skips/acceptances
      test <- res[, .(hitmix = length(unique(hit)) == 1), by = hit_grp]
      stopifnot(all(test$hitmix))
    } else {
      res[, fin_grp := NA]
      res[, hit_grp := NA]
    }
  }
  # return fewer columns to save memory
  returncols <- names(res)
  if (only_hits) {
    res <- droplevels(res[(hit), .SD, .SDcols = returncols])
  }
  return(res[, .SD, .SDcols = returncols])
}

#' Make a plot of the nodes
#'
#' Given the results of the splitting and testing algorithm in the form of a
#' graph from [make_results_tree], make a node level data set for use in
#' reporting results in terms of a binary tree graph. This does not print or
#' plot the graph. You'll need to do that with the resulting object.
#'
#' @param res_graph A tidygraph object produced from make_results_tree

#' @param remove_na_p A logical indicating whether the graph should include
#' nodes/leaves that were not tested. Default (TRUE) is to remove them. When
#' remove_na_p is FALSE, the graph might look strange since some blocks will
#' not have a known position in the graph  (the graph is fixed, but not
#' specified by the find_blocks function when a node or block is not visited
#' for testing.)

#' @return A ggraph object
#' @import ggraph
#' @import ggplot2
#' @export
make_results_ggraph <- function(res_graph, remove_na_p = TRUE) {
  if (remove_na_p) {
    res_graph <- res_graph %>%
      activate(nodes) %>%
      filter(!is.na(p))
  }
  res_g <- ggraph(res_graph, layout = "tree") +
    geom_edge_diagonal(colour = "black") +
    geom_node_label(aes(label = label, colour = hit),
      repel = FALSE, show.legend = FALSE, label.r = unit(0.5, "lines"),
      label.padding = unit(.01, "lines"), label.size = 0
    ) +
    theme(legend.position = "none") +
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  return(res_g)
}

#' Make a node-level dataset from a block-level dataset
#'
#' Given the results of the splitting and testing algorithm, make a node level
#' data set for use in reporting results and as input to ggraph for
#' visualization in terms of a tree graph.
#'
#' @param orig_res data.table from find_blocks(); must include elements such as
#'   - biggrp    (dot-sep lineage, may be truncated),
#'   - p1,p2,… and a pfinal*,
#'   - alpha1, alpha2, …
#' @param block_id   name of your block ID column (e.g. "bF")
#' @param node_label optional name of a descriptive label column
#' @param return_what a character vector containing "all", "graph" (a tbl_graph
#' object with nodes and edges), "nodes" (a data.table with node level
#' information), "test_summary" (a data.table object with one row indicating
#' false and true discoveries, etc.)

#' @param truevar_name the optional name of a column recording the true
#' treatment effect (used here to find blocks where the true effect is 0 or
#' not). In some simulations we have a column called nonnull which is TRUE if
#' that block or node has a non-zero effect and FALSE if the block or node has
#' a truly zero effect. So, truevar_name can be "nonnull"

#' @return a list that can contain nodes, a tbl_graph object, and/or a test_summary
#' @importFrom stringi stri_split_fixed stri_sub
#' @importFrom tidygraph tbl_graph centrality_degree node_is_adjacent activate
#' @import tidygraph data.table
#' @export
make_results_tree <- function(orig_res, block_id, node_label = NULL, return_what = "all", truevar_name = NULL) {
  dt <- copy(orig_res)

  # assign a unique leaf ID per input row
  stopifnot(all.equal(length(unique(dt[[block_id]])), nrow(dt)))
  dt[, leaf_id := as.character(get(block_id))]
  ## hashing because some block_ids use numbers like "1" which is also the number of the root node
  dt[leaf_id == "1", leaf_id := nodeidfn(leaf_id)]

  # split the biggrp chains
  ancestry <- stri_split_fixed(dt$biggrp, ".", simplify = TRUE)
  ancestry[ancestry == ""] <- NA_character_
  max_depth <- ncol(ancestry)

  node_lvls <- paste0("node", seq_len(max_depth))
  dt[, (node_lvls) := as.data.table(ancestry)]
  ## For any missing leaves, fill in with the leaf_id for the leaf level
  dt[, (node_lvls[max_depth]) := ifelse(is.na(get(node_lvls[max_depth])), leaf_id, get(node_lvls[max_depth]))]

  # stopifnot(all.equal(length(unique(dt[[node_lvls[max_depth]]])), nrow(dt)))

  # detect p-columns: p1,p2,… then your pfinal* as the last one
  p_normals <- grep("^p[0-9]+$", names(dt), value = TRUE)
  p_final <- grep("^pfinal", names(dt), value = TRUE)
  if (length(p_final) == 0) stop("No pfinal* column found.")
  pcols <- c(p_normals, p_final[1])

  # detect alpha-columns similarly
  a_cols <- grep("^alpha[0-9]+$", names(dt), value = TRUE)

  if (!is.null(truevar_name)) {
    dt[, nonnull := get(truevar_name) != 0]
  } else {
    dt[, nonnull := NA]
  }

  if (is.null(node_label)) {
    dt[, node_label := NA_character_]
  } else {
    dt[, node_label := get(node_label)]
  }

  # Convert from what is in essence wide-form to long-form (which is node-form)
  node_list <- vector("list", max_depth)

  for (d in seq_len(max_depth)) {
    this_pcol <- pcols[d]
    this_acol <- a_cols[d]
    ids <- dt[[node_lvls[d]]]
    if (any(is.na(ids))) {
      ## If this pasting becomes unwieldy we can use nodeidfn() to replace it with a hash
      dt[is.na(get(node_lvls[d])), (node_lvls[d]) := paste(unique(leaf_id), collapse = ","), by = eval(node_lvls[d - 1])]
    }

    res <- dt[, .(
      depth = d,
      parent_name = unique(get(node_lvls[max(1, d - 1)])),
      p = unique(get(this_pcol)),
      a = unique(get(this_acol)),
      ## Remember that any ancester of a non-null leaf is non-null
      nonnull = any(nonnull, na.rm = TRUE),
      blocks = unique(paste(unique(leaf_id), collapse = ",")),
      node_label = unique(paste(unique(node_label), collapse = ",")),
      num_leaves = .N
    ),
    by = .(name = ids)
    ]

    ## Assign node numbers: we cannot guarantee a complete k-ary tree so we use
    ## the results of the previous loop. We mostly have these for labeling the graph. We might do something else.
    ## Since these here do not follow the parent child relationships that we'd like.
    if (d == 1) {
      res[, node_number := 1]
      res[, parent_name := NA]
    } else {
      offset <- max(node_list[[d - 1]]$node_number)
      k <- nrow(res)
      res[, node_number := offset + seq_len(k)]
    }

    node_list[[d]] <- res
  }
  nodes_dt <- rbindlist(node_list, use.names = TRUE)

  ## Can't calculate errors and discoveries without knowing the truth aka
  ## having a variable that records whether a node/block is not null or not.

  if (any(!is.na(nodes_dt$nonnull))) {
    # Record testing results
    num_nodes <- nrow(nodes_dt)
    # the block level dataset has one row for each leaf or block
    num_leaves <- nrow(dt)
    num_nodes_tested <- sum(!is.na(nodes_dt$p))
    num_nonnull_nodes_tested <- sum(!is.na(nodes_dt$p) & nodes_dt$nonnull)
    ## A discovery is a rejection
    node_rejections <- nodes_dt[!is.na(p), sum(p <= a, na.rm = TRUE)]
    ## Record rejections info
    node_any_false_rejection <- nodes_dt[nonnull == FALSE & !is.na(p), any(p <= a)]
    node_false_rejection_prop <- nodes_dt[nonnull == FALSE & !is.na(p), mean(p <= a)]
    node_num_false_rejections <- nodes_dt[nonnull == FALSE & !is.na(p), sum(p <= a)]
    ## Now false discovery prop (element in FDR calc). Denominator is rejections.
    node_false_discovery_prop <- nodes_dt[nonnull == FALSE & !is.na(p), sum(p <= a) / max(1, node_rejections)]

    ## prop of nodes (including leaves) where p should be less than or equal to a
    ## prop of true discoveries among the possible true discoveries
    node_true_discoveries <- nodes_dt[nonnull == TRUE & !is.na(p), sum(p <= a)]
    node_power <- nodes_dt[nonnull == TRUE & !is.na(p), mean(p <= a)]

    ## If max_depth=1, then we are only testing the root node
    ## And we are excluding for now the idea that the root is the leaf --- that is just the case of a single test.
    num_leaves_tested <- (max_depth > 1) * sum(nodes_dt[num_leaves == 1, !is.na(p)])
    num_nonnull_leaves_tested <- (max_depth > 1) * sum(nodes_dt[num_leaves == 1 & nonnull == TRUE, !is.na(p)])
    if (num_leaves_tested > 0) {
      leaf_power <- nodes_dt[num_leaves == 1 & nonnull == TRUE & !is.na(p), mean(p <= a, na.rm = TRUE)]
      leaf_rejections <- nodes_dt[num_leaves == 1 & !is.na(p), sum(p <= a, na.rm = TRUE)]
      leaf_true_discoveries <- nodes_dt[num_leaves == 1 & nonnull == TRUE & !is.na(p), sum(p <= a, na.rm = TRUE)]
      leaf_any_false_rejection <- nodes_dt[num_leaves == 1 & nonnull == FALSE & !is.na(p), any(p <= a)]
      leaf_false_rejection_prop <- nodes_dt[num_leaves == 1 & nonnull == FALSE & !is.na(p), mean(p <= a, na.rm = TRUE)]
      leaf_false_discovery_prop <- nodes_dt[num_leaves == 1 & nonnull == FALSE & !is.na(p), sum(p <= a, na.rm = TRUE) / max(1, leaf_rejections)]
    } else {
      leaf_power <- 0
      leaf_rejections <- 0
      leaf_true_discoveries <- 0
      leaf_any_false_rejection <- 0
      leaf_false_rejection_prop <- 0
      leaf_false_discovery_prop <- 0
    }

    test_summary <- data.table(
      num_nodes = num_nodes,
      num_leaves = num_leaves,
      num_nodes_tested = num_nodes_tested,
      num_nonnull_nodes_tested = num_nonnull_nodes_tested,
      node_rejections = node_rejections,
      node_any_false_rejection = node_any_false_rejection,
      node_false_rejection_prop = node_false_rejection_prop,
      node_num_false_rejections = node_num_false_rejections,
      node_false_discovery_prop = node_false_discovery_prop,
      node_true_discoveries = node_true_discoveries,
      node_power = node_power,
      num_leaves_tested = num_leaves_tested,
      num_nonnull_leaves_tested = num_nonnull_leaves_tested,
      leaf_power = leaf_power,
      leaf_rejections = leaf_rejections,
      leaf_true_discoveries = leaf_true_discoveries,
      leaf_any_false_rejection = leaf_any_false_rejection,
      leaf_false_rejection_prop = leaf_false_rejection_prop,
      leaf_false_discovery_prop = leaf_false_discovery_prop
    )
  } else {
    test_summary <- NA
  }

  nodes_dt[, hit := p <= a]

  if (any(c("all", "graph") %in% return_what)) {
    ### set up the tidygraph style object
    edges_dt <- nodes_dt[!is.na(parent_name), .(name, parent_name)]
    setnames(edges_dt, c("parent_name", "name"), c("from", "to"))
    res_graph <- tbl_graph(nodes = nodes_dt, edges = edges_dt)

    ## Abbreviate the block name string variable
    res_graph <- res_graph %>%
      activate(nodes) %>%
      mutate(shortbf = ifelse(nchar(blocks) > 6, paste0(stri_sub(blocks, 1, 5), "..."), blocks))
    ## Make names etc for ease in plotting later
    if (length(unique(nodes_dt$a)) > 1) {
      res_graph <- res_graph %>%
        activate(nodes) %>%
        mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf,
          "\n # Blocks=", num_leaves, "\n p=", round(p, 3), ",a=", round(a, 3),
          sep = ""
        ))
    } else {
      res_graph <- res_graph %>%
        activate(nodes) %>%
        mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf,
          "\n # Blocks=", num_leaves, "\n p=", round(p, 3),
          sep = ""
        ))
    }
  } else {
    res_graph <- NULL
  }

  results <- list(nodes = nodes_dt, graph = res_graph, test_summary = test_summary)
  if ("all" %in% return_what) {
    return_what <- c("nodes", "graph", "test_summary")
  }
  return(results[return_what])
}

### DEPRECATING THE BELOW
### #' Make a node level tree object of the results of nested testing
### #'
### #' Given the results of the splitting and testing algorithm, make a node level
### #' data set for use in reporting results and as input to ggraph for
### #' visualization in terms of a tree graph.
### #'
### #' @param orig_res a results data.table output from the \code{\link{find_blocks}} function.
### #' @param blockid is a character name for the variable containing the block id information
### #' @param node_label is a character name for a variable containing a descriptive label for the blocks.
### #' @return A tbl_graph and igraph object with nodes and edges
### #' @importFrom stringi stri_split_fixed stri_sub
### #' @importFrom tidygraph tbl_graph centrality_degree node_is_adjacent activate
### #' @import tidygraph
### #' @importFrom data.table melt
### #' @export
### make_results_tree <- function(orig_res, blockid = "bF", node_label = NULL) {
###   # We have to make a node level data set and an edge level data set in order to define the graph
###   res <- copy(orig_res)
###   ## for testins
###   # res <- res_fwer
###   # blockid <- "bF"
###   res_nodeids <- stri_split_fixed(as.character(res$biggrp), pattern = ".", simplify = TRUE)
###   res_nodeids[res_nodeids == ""] <- NA
###   # class(res_nodeids) <- "integer"
###   res[, paste("nodenum", 1:ncol(res_nodeids), sep = "") := lapply(1:ncol(res_nodeids), function(i) {
###     res_nodeids[, i]
###   })]
###   pnms <- sort(grep("^p[0-9]", names(res), value = TRUE))
###   anms <- sort(grep("^alpha[0-9]", names(res), value = TRUE))
###   nodenums <- sort(grep("^nodenum[0-9]", names(res), value = TRUE))
###   # TODO Right now we go from wide to long to node. It would be nicer to go directly from wide to node.
###
###   ## Not all results will have either a node_label defined or a nonnull vector, so we create them with NULL if needed
###   if (is.null(node_label)) {
###     res[, node_label := "NULL"]
###   } else {
###     res[, node_label := get(node_label)]
###   }
###   if (!any(names(res) == "nonnull")) {
###     res[, nonnull := NULL]
###   }
###   longnms <- c("biggrp", blockid, "nodenum_current", "nodenum_prev", "node_label", "nonnull")
###   reslong <- melt(res,
###     id = longnms,
###     measure.vars = list(p = pnms, a = anms, nodenum = nodenums),
###     variable.name = "depth"
###   )
###   reslong$depth <- as.numeric(as.character(reslong$depth))
###   reslong$bFC <- as.character(reslong[[blockid]])
###   reslong <- droplevels(reslong[!is.na(nodenum) & !is.na(p), ])
###   ## Now collapse down to the node level
###   res_nodes_df <- reslong[, .(
###     p = unique(p),
###     a = unique(a),
###     bF = paste(as.character(unlist(sort(get(blockid)))), collapse = ","),
###     depth = unique(depth),
###     # carry through any signal from leaves
###     node_nonnull = any(nonnull, na.rm = TRUE),
###     node_label = paste(as.character(unlist(sort(node_label))), collapse = ",")
###   ), by = nodenum]
###   res_nodes_df$name <- res_nodes_df$nodenum
###   res_nodes_df$num_blocks <- stri_count_fixed(res_nodes_df$bF, ",") + 1
###
###   # Make an edge data.frame
###   res_edges_lst <- list()
###   if (nrow(res_nodes_df) == 1) {
###     res_edges_lst[[1]] <- data.table(from = "1", to = "1")
###   } else {
###     for (i in 1:(ncol(res_nodeids) - 1)) {
###       res_edges_lst[[i]] <- data.table(from = res_nodeids[, i], to = res_nodeids[, i + 1])
###     }
###   }
###   res_edges_df <- unique(na.omit(rbindlist(res_edges_lst)))
###   res_nodes_df <- merge(res_nodes_df, res_edges_df,
###     by.x = "nodenum", by.y = "to",
###     sort = FALSE, all.x = TRUE
###   )
###   setnames(res_nodes_df, "from", "parent_name")
###
###   # Now define the graph using the node data set and the edges dataset.
###   res_graph <- tbl_graph(nodes = res_nodes_df, edges = res_edges_df)
###
###   # And use the graph relations to calculate whether a test at a given place in the tree is a discovery or not
###
###   # first way to detect is leaf with p=<alpha and second way as parent of all non-sig leaves
###   # leaf is a single experimental block here at the end of the tree. A node
###   # that consists of a single block.
###
###   res_graph <- res_graph %>%
###     activate(nodes) %>%
###     mutate(
###       out_degree = centrality_degree(mode = "out"),
###       is_leaf = node_is_leaf(),
###       is_leaf_single_block = (out_degree == 0 & depth > 1 & num_blocks == 1)
###     ) %>%
###     group_by(parent_name) %>%
###     mutate(leaf_child_all_not_sig = all(p > a & is_leaf)) %>%
###     ungroup()
###
###   res_graph <- res_graph %>%
###     activate(nodes) %>%
###     mutate(
###       leaf_hit = (p <= a & is_leaf_single_block),
###       is_leaf_parent = node_is_adjacent(to = is_leaf_single_block, mode = "in", include_to = FALSE),
###       is_leaf_parent2 = name %in% unique(parent_name[is_leaf_single_block]),
###       num_desc = local_size(order = graph_order(), mode = "out", mindist = 1),
###       is_cut = node_is_cut(),
###       parent_of_all_notsig_leaves = node_is_adjacent(to = leaf_child_all_not_sig, mode = "in", include_to = FALSE)
###     )
###   stopifnot(all.equal(res_graph$is_leaf_parent, res_graph$is_leaf_parent2))
###
###   ## ## the is_cut nodes are those at the base of the tree --- no further splitting
###   ## ## some of them are leaves (individual blocks) and others are groups of blocks (not leaves)
###
###   ## res_graph %>%
###   ##   filter(!is_leaf_single_block) %>%
###   ##   select(nodenum, depth, num_blocks, parent_name, out_degree, is_leaf, is_leaf_parent, is_leaf_parent2, p, num_desc, dist_to_leaf, is_cut) %>%
###   ##   as_tibble() %>%
###   ##   print(n = 100)
###
###   ## We use group_hit for indirect discovery (i.e. we can reject the null of no effects in any of these blocks, but not in one or the other block
###   res_graph <- res_graph %>%
###     activate(nodes) %>%
###     mutate(
###       group_hit = (p <= a & parent_of_all_notsig_leaves),
###       hit = group_hit | leaf_hit
###     )
###
###   ## Abbreviate the block name string variable
###   res_graph <- res_graph %>%
###     activate(nodes) %>%
###     mutate(shortbf = ifelse(nchar(bF) > 6, paste0(stri_sub(bF, 1, 5), "..."), bF))
###   ## Make names etc for ease in graphing later
###   if (length(unique(res_nodes_df$a)) > 1) {
###     res_graph <- res_graph %>%
###       activate(nodes) %>%
###       mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf, "\n # Blocks=", num_blocks, "\n p=", round(p, 3), ",a=", round(a, 3), sep = ""))
###   } else {
###     res_graph <- res_graph %>%
###       activate(nodes) %>%
###       mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf, "\n # Blocks=", num_blocks, "\n p=", round(p, 3), sep = ""))
###   }
###
###   return(res_graph)
### }
