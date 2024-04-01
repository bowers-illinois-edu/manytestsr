# Functions for organizing and graphing results

#' Return detected blocks plus info
#'
#' Given the results of the splitting and testing algorithm, report on the blocks
#'  where the null of no effects could be rejected at level alpha. Currently calculates rejections using an FWER style criteria (p of a node = max of all previous nodes) if the final alphas are all the same as the scalar alpha OR if fwer=TRUE.
#' @param orig_res results data.table output from the \code{\link{findBlocks}} function.
#' @param fwer (default is TRUE) means that a block is detected (or not) using the maximum p-value associated with the
#' block (or the groups containing that block). fwer=FALSE to detect blocks (or groups of blocks) using FDR control.
#' @param alpha Is the false positive rate used for detecting an effect if it is constant (i.e. not an FDR-style approach).
#' @param only_hits (default FALSE) returns only the detected blocks instead of all of them
#' @param blockid Name of block variable (the blocking variable is a factor)
#' @return A data.table adding a column \code{hit} to the \code{res} data.table indicating a "hit" or detection for that block (or group of blocks)
#' @importFrom stringi stri_count_fixed stri_split_fixed
#' @import data.table
#' @export
report_detections <- function(orig_res, fwer = TRUE, alpha = .05, only_hits = FALSE, blockid = "blockF") {
  res <- copy(orig_res)
  # For the output from adjust_block_tests (the bottom-up or test all blocks method)
  if (length(grep("biggrp", names(res))) == 0) {
    # This is for the bottom-up/test every block method
    # max_p are the adjusted p-values so we can use alpha=.05 for error rate control (hoping it is fdr and not fwer)
    res[, hit := max_p <= alpha]
    res[, hit_grp := nodenum_current]
    if (all(!res$hit)) {
      res[, hit_grp := NA]
    }
  } else {
    # For the splitting based methods (output  from findBlocks)
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
    # If the block is a terminal node and p<=alpha then this is a hit or detected effect
    # If the block is a terminal node and p>alpha for both it and all other terminal nodes
    #  then we have a hit or detected effect for *all* terminal nodes but cannot distinguish between them.
    # I think max_p should be p at maxdepth for the FDR algorithm and maximum p overall (or pfinalb) for the FWER algo.
    # could have done this first     res[, max_alpha := get(paste0("alpha", maxdepth)), by = seq_len(nrow(res))]
    # But I think the following is faster
    if (fwer | all(res$parent_alpha == alpha)) {
      res[, max_p := pfinalb]
      res[, max_alpha := alpha]
    } else {
      res[, max_p := get(paste0("p", maxdepth)), by = seq_len(nrow(res))]
      res[, max_alpha := get(paste0("alpha", maxdepth)), by = seq_len(nrow(res))]
    }
    # A detection on a single block is scored if the final p <= alpha for a node containing a single block (i.e. a leaf)
    res[, single_hit := max_p <= max_alpha & blocksbygroup == 1]

    # A detection is also scored if all leaves have p > alpha but the parent
    # has p <= alpha: this is a grouped detection with multiple blocks.

    res[, group_hit := fifelse(!single_hit & (all(max_p > max_alpha) & (fin_parent_p <= parent_alpha) & blocksbygroup >= 1), TRUE, FALSE), by = fin_parent]

    # Also a group hit can be scored (an effect detected within a group of
    # blocks) if there are multiple blocks in a final node and that test is p<a

    res[, group_hit2 := fifelse(blocksbygroup > 1 & all(max_p <= max_alpha), TRUE, FALSE), by = fin_nodenum]
    res[, hit := single_hit | group_hit | group_hit2]
    if (any(res$hit)) {
      res[(hit), fin_grp := fifelse(single_hit | group_hit2, fin_nodenum, fin_parent)]
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
  # Later return fewer columns to save memory
  returncols <- names(res)
  if (only_hits) {
    res <- droplevels(res[(hit), .SD, .SDcols = returncols])
  }
  return(res[, .SD, .SDcols = returncols])
}


#' Make a node level tree object of the results of nested testing
#'
#' Given the results of the splitting and testing algorithm, make a node level
#' data set for use in reporting results and as input to ggraph for
#' visualization in terms of a tree graph.
#'
#' @param orig_res results data.table output from the \code{\link{findBlocks}} function.
#' @param blockid is Is a character name for the variable containing the block id information
#' @return A tbl_graph and igraph object with nodes and edges
#' @importFrom stringi stri_split_fixed stri_sub
#' @importFrom tidygraph tbl_graph centrality_degree node_is_adjacent activate
#' @import tidygraph
#' @importFrom data.table melt
#' @export
make_results_tree <- function(orig_res, blockid = "bF") {
  # We have to make a node level data set and an edge level data set in order to define the graph
  res <- copy(orig_res)
  ## for testins
  # res <- res_fwer
  # blockid <- "bF"
  res_nodeids <- stri_split_fixed(as.character(res$biggrp), pattern = ".", simplify = TRUE)
  res_nodeids[res_nodeids == ""] <- NA
  # class(res_nodeids) <- "integer"
  res[, paste("nodenum", 1:ncol(res_nodeids), sep = "") := lapply(1:ncol(res_nodeids), function(i) {
    res_nodeids[, i]
  })]
  pnms <- sort(grep("^p[0-9]", names(res), value = TRUE))
  anms <- sort(grep("^alpha[0-9]", names(res), value = TRUE))
  nodenums <- sort(grep("^nodenum[0-9]", names(res), value = TRUE))
  # Right now we go from wide to long to node. It would be nicer to go directly from wide to node.
  reslong <- melt(res,
    id = c("biggrp", blockid, "nodenum_current", "nodenum_prev"),
    measure.vars = list(p = pnms, a = anms, nodenum = nodenums),
    variable.name = "depth"
  )
  reslong$depth <- as.numeric(as.character(reslong$depth))
  reslong$bFC <- as.character(reslong[[blockid]])
  reslong <- droplevels(reslong[!is.na(nodenum) & !is.na(p), ])

  res_nodes_df <- reslong[, .(
    p = unique(p),
    a = unique(a),
    bF = paste(as.character(unlist(sort(get(blockid)))), collapse = ","),
    depth = unique(depth)
  ), by = nodenum]
  res_nodes_df$name <- res_nodes_df$nodenum
  res_nodes_df$num_blocks <- stri_count_fixed(res_nodes_df$bF, ",") + 1

  # Make an edge data.frame
  res_edges_lst <- list()
  if (nrow(res_nodes_df) == 1) {
    res_edges_lst[[1]] <- data.table(from = "1", to = "1")
  } else {
    for (i in 1:(ncol(res_nodeids) - 1)) {
      res_edges_lst[[i]] <- data.table(from = res_nodeids[, i], to = res_nodeids[, i + 1])
    }
  }
  res_edges_df <- unique(na.omit(rbindlist(res_edges_lst)))
  res_nodes_df <- merge(res_nodes_df, res_edges_df,
    by.x = "nodenum", by.y = "to",
    sort = FALSE, all.x = TRUE
  )
  setnames(res_nodes_df, "from", "parent_name")

  # Now define the graph using the node data set and the edges dataset.
  res_graph <- tbl_graph(nodes = res_nodes_df, edges = res_edges_df)

  # first way to detect is leaf with p=<alpha and second way as parent of all non-sig leaves
  # leaf is a single experimental block here at the end of the tree. A node
  # that consists of a single block.

  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(
      out_degree = centrality_degree(mode = "out"),
      is_leaf = node_is_leaf(),
      is_leaf_single_block = (out_degree == 0 & depth > 1 & num_blocks == 1)
    ) %>%
    group_by(parent_name) %>%
    mutate(leaf_child_all_not_sig = all(p > a & is_leaf)) %>%
    ungroup()

  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(
      leaf_hit = (p <= a & is_leaf_single_block),
      is_leaf_parent = node_is_adjacent(to = is_leaf_single_block, mode = "in", include_to = FALSE),
      is_leaf_parent2 = name %in% unique(parent_name[is_leaf_single_block]),
      num_desc = local_size(order = graph_order(), mode = "out", mindist = 1),
      is_cut = node_is_cut(),
      parent_of_all_notsig_leaves = node_is_adjacent(to = leaf_child_all_not_sig, mode = "in", include_to = FALSE)
    )
  stopifnot(all.equal(res_graph$is_leaf_parent, res_graph$is_leaf_parent2))

  ## single_block_leaf_names <- res_graph %>%
  ##   activate(nodes) %>%
  ##   filter(is_leaf_single_block) %>%
  ##   pull(nodenum)

  ## distances(graph = as.igraph(res_graph), v = res_graph %>% activate(nodes) %>% filter(is_leaf_single_block), to = res_graph %>% activate(nodes) %>% filter(!is_leaf_single_block))

  ## all_dists <- distances(as.igraph(res_graph), mode = "out")
  ## tmp <- all_dists[!(row.names(all_dists) %in% single_block_leaf_names), single_block_leaves_names]
  ## tmp_max_dist <- apply(tmp, 1, function(x) {
  ##   newx <- x[!is.infinite(x)]
  ##   max(newx, na.rm = TRUE)
  ## })

  ## ## From https://stackoverflow.com/questions/69496134/how-to-get-all-leaf-nodes-from-a-directed-subtree-using-igraph-in-r
  ## library(igraph)
  ## f <- function(g, r) {
  ##   names(V(g))[is.finite(distances(g, r, mode = "out")) & degree(g) == 1]
  ## }

  ## fun <- function(graph, node) {
  ##   path <- ego(graph, order = length(V(graph)), nodes = node, mode = "out")
  ##   nms <- names(path[[1]])
  ##   ## nms[ego_size(graph, order=1, nodes=nms, mode="out", mindist=1) == 0]
  ##   nms[degree(graph, v = nms, mode = "out") == 0]
  ## }

  ## res_igraph <- as.igraph(res_graph)
  ## fun(res_graph, node = single_block_leaf_names[1])
  ## f(res_graph, single_block_leaves_names[1])

  ## distances(g = res_igraph, v = V(single_block_leaf_names), to = V(res_igraph))

  ## ## the is_cut nodes are those at the base of the tree --- no further splitting
  ## ## some of them are leaves (individual blocks) and others are groups of blocks (not leaves)

  ## res_graph %>%
  ##   filter(!is_leaf_single_block) %>%
  ##   select(nodenum, depth, num_blocks, parent_name, out_degree, is_leaf, is_leaf_parent, is_leaf_parent2, p, num_desc, dist_to_leaf, is_cut) %>%
  ##   as_tibble() %>%
  ##   print(n = 100)

  ## We use group_hit for indirect discovery (i.e. we can reject the null of no effects in any of these blocks, but not in one or the other block
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(
      group_hit = (p <= a & parent_of_all_notsig_leaves),
      hit = group_hit | leaf_hit
    )

  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(shortbf = ifelse(nchar(bF) > 6, paste0(stri_sub(bF, 1, 5), "..."), bF))
  if (length(unique(res_nodes_df$a)) > 1) {
    res_graph <- res_graph %>%
      activate(nodes) %>%
      mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf, "\n # Blocks=", num_blocks, "\n p=", round(p, 3), ",a=", round(a, 3), sep = ""))
  } else {
    res_graph <- res_graph %>%
      activate(nodes) %>%
      mutate(label = paste("Node:", stri_sub(name, 1, 4), "\n Name:", shortbf, "\n # Blocks=", num_blocks, "\n p=", round(p, 3), sep = ""))
  }

  return(res_graph)
}

#' Make a plot of the nodes
#'
#' Given the results of the splitting and testing algorithm in the form of a
#' graph from [make_results_tree], make a node level data set for use in reporting results in terms of a binary tree graph. This does not print or plot the graph. You'll need to do that with the resulting object.
#' @param res_graph A tidygraph object produced from make_results_tree
#' @return A ggraph object
#' @import ggraph
#' @import ggplot2
#' @export
make_results_ggraph <- function(res_graph) {
  # require(ggraph)
  # require(tidygraph)
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
