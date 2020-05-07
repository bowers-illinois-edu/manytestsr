## Functions for organizing and graphing results


##' Return detected blocks plus info
##'
##' Given the results of the splitting and testing algorithm, report on the blocks
##'  where the null of no effects could be rejected at level alpha.
##' @param orig_res results data.table output from the \code{\link{findBlocks}} function.
##' @param fwer (default is TRUE) means that a block is detected (or not) using the maximum p-value associated with the
##' block (or the groups containing that block). fwer=FALSE to detect blocks (or groups of blocks) using FDR control.
##' @param only_hits (default FALSE) returns only the detected blocks
##' @param autofwer If fwer=TRUE but alpha varies, return the fdr based report.
##' @return A data.table adding a column \code{hit} to the \code{res} data.table indicating a "hit" or detection for that block (or group of blocks)
##' @importFrom stringi stri_count_fixed stri_split_fixed
##' @export
report_detections <- function(orig_res, fwer = TRUE, alpha = .05, only_hits = FALSE, autofwer=TRUE) {
  require(stringi) ## comment out for production
  res <- copy(orig_res)
  res_nodeids <- stri_split_fixed(as.character(res$biggrp), pattern = ".", simplify = TRUE)
  class(res_nodeids) <- "integer"
  ## Extract final node number (findBlocks uses the canonical binary tree node numbering system)
  res_fin_nodenum <- apply(res_nodeids, 1, function(x) {
    max(x, na.rm = TRUE)
  })
  res[, fin_nodenum := res_fin_nodenum]
  res[, fin_parent := floor(fin_nodenum / 2)]
  ## Maximum tree depth for a node encoded in the biggrp string: basically number of dots+1 or number of node numbers
  ## (where node numbers are separated by a dot)
  res[, maxdepth := stri_count_fixed(biggrp, ".") + 1]
  if (all(res$maxdepth == 1)) {
    ## When the algorithm stops at the first level there are no parents. There is just the root node. This is an attempt at a work around
    res[, fin_parent_p := p1]
    res[, parent_alpha := alpha1]
  } else {
    res[, fin_parent_p := get(paste0("p", maxdepth - 1)), by = seq_len(nrow(res))]
    res[, parent_alpha := get(paste0("alpha", maxdepth - 1)), by = seq_len(nrow(res))]
  }
  ## If the block is a terminal node and p<alpha then this is a hit or detected effect
  ## If the block is a terminal node and p>alpha for both it and the other of the two terminal nodes then we have a hit or detected effect for *both* nodes but cannot distinguish between them.
  ## I think max_p should be p at maxdepth for the FDR algorithmn and maximum p overall (or pfinalb) for the FWER algo.
  ## could have done this first     res[, max_alpha := get(paste0("alpha", maxdepth)), by = seq_len(nrow(res))]
  ## But I think the following is faster
  if (fwer & !autofwer) {
    res[, max_p := pfinalb]
    res[, max_alpha := alpha]
    # res[, parent_alpha := alpha]
  } else {
    res[, max_p := get(paste0("p", maxdepth)), by = seq_len(nrow(res))]
    res[, max_alpha := get(paste0("alpha", maxdepth)), by = seq_len(nrow(res))]
  }
  ## A detection is scored if p < alpha for a node containing a single block (i.e. a leaf)
  res[, single_hit := max_p < max_alpha & blocksbygroup == 1]
  ## A detection is also scored if the two leaves have p > alpha but the parent has p < alpha: this is a grouped detection with two blocks.
  res[(!single_hit), group_hit := all(max_p > max_alpha) & fin_parent_p < parent_alpha & length(unique(fin_nodenum)) == 2, by = fin_parent]
  res[, hit := single_hit | group_hit]
  ## This records the group of blocks or individual blocks detected
  if (any(res$hit)) {
    res[(hit), fin_grp := ifelse(single_hit, fin_nodenum, fin_parent)]
    res[(hit), hit_grp := factor(fin_grp)] #, labels = 1:length(unique(fin_grp)))]
    res[!(hit), hit_grp := factor(nodenum_current)]
    ## Make sure that no  hit_grp includes *both* detections and misses/skips/acceptances
    test <- res[,.(hitmix= length(unique(hit))==1),by=hit_grp]
    stopifnot(all(test$hitmix))
  } else {
    res[, fin_grp := NA]
    res[, hit_grp := NA]
  }
  returncols <- names(res) ##c("bF", "biggrp", "hit", "hit_grp", "fin_grp", "max_p", "max_alpha", "fin_parent_p", "parent_alpha", grep("^ate", names(res), value = TRUE))
  if (only_hits) {
    res <- droplevels(res[(hit), .SD, .SDcols = returncols])
  }
  return(res[, .SD, .SDcols = returncols])
}


##' Make a node level binary tree object
##'
##' Given the results of the splitting and testing algorithm, make a node level data set for use in reporting results in terms of a binary tree graph.
##' @param orig_res results data.table output from the \code{\link{findBlocks}} function.
##' @param fwer Means that alpha is fixed otherwise use the alpha that is calculated at each step of the splitting procedure
##' @return A tbl_graph and igraph object with nodes and edges
##' @importFrom stringi stri_split_fixed
##' @importFrom tidygraph tbl_graph centrality_degree node_is_adjacent
##' @importFrom data.table melt
##' @export
make_tree <- function(orig_res) {
  require(ggraph)
  require(tidygraph)
  ## We have to make a node level data set and an edge level data set in order to define the graph
  res <- copy(orig_res)
  res_nodeids <- stri_split_fixed(as.character(res$biggrp), pattern = ".", simplify = TRUE)
  class(res_nodeids) <- "integer"
  res[, paste("nodenum", 1:ncol(res_nodeids), sep = "") := lapply(1:ncol(res_nodeids), function(i) {
    res_nodeids[, i]
  })]
  pnms <- sort(grep("^p[0-9]", names(res), value = TRUE))
  anms <- sort(grep("^alpha[0-9]", names(res), value = TRUE))
  nodenums <- sort(grep("^nodenum[0-9]", names(res), value = TRUE))
  ## Right now we go from wide to long to node. It would be nicer to go directly from wide to node.
  reslong <- melt(res,
    id = c("biggrp", "bF"),
    measure.vars = list(p = pnms, a = anms, nodenum = nodenums),
    variable.name = "depth"
  )
  reslong$depth <- as.numeric(as.character(reslong$depth))
  reslong$bFC <- as.character(reslong$bF)
  reslong <- droplevels(reslong[!is.na(nodenum) & !is.na(p), ])
  res_nodes_df <- reslong[, .(p = unique(p), a = unique(a), bF = paste(as.character(unlist(sort(bF))), collapse = ","), depth = unique(depth)), by = nodenum]
  res_nodes_df$name <- res_nodes_df$nodenum
  ## Make an edge data.frame
  res_edges_lst <- list()
  for (i in 1:(ncol(res_nodeids) - 1)) {
    res_edges_lst[[i]] <- data.table(from = res_nodeids[, i], to = res_nodeids[, i + 1])
  }
  res_edges_df <- unique(na.omit(rbindlist(res_edges_lst)))
  res_edges_df[, c("to", "from") := .(as.character(to), as.character(from))]
  # res_edges_df <- res_edges_df %>% lazy_dt() %>% mutate_at(vars(to, from), as.character) %>% as.data.table()
  ## Now define the graph using the node data set and the edges dataset.
  res_graph <- tbl_graph(nodes = res_nodes_df, edges = res_edges_df)
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(out_degree = centrality_degree(mode = "out"))
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(parent_name = as.character(floor(as.numeric(name) / 2)))
  res_graph <- res_graph %>%
    activate(nodes) %>%
    group_by(parent_name) %>%
    mutate(both_ns = all(p > a)) %>%
    ungroup()
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(parent_of_ns = node_is_adjacent(to = both_ns, include_to = FALSE, mode = "in"))
  ## first way to detect is leaf with p<alpha and second way as parent of two non-sig leaves
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(hit = (out_degree == 0 & p <= a) | parent_of_ns)
  res_graph <- res_graph %>%
    activate(nodes) %>%
    mutate(shortbf = ifelse(nchar(bF) > 6, paste0(substr(bF, 1, 5), "..."), bF))
  if (length(unique(res_nodes_df$a)) > 1) {
    res_graph <- res_graph %>%
      activate(nodes) %>%
      mutate(label = paste("Node:", name, "\n Blocks:", shortbf, "\n p=", round(p, 3), ",a=", round(a, 3), sep = ""))
  } else {
    res_graph <- res_graph %>%
      activate(nodes) %>%
      mutate(label = paste("Node:", name, "\n Blocks:", shortbf, "\n p=", round(p, 3), sep = ""))
  }
  return(res_graph)
}

##' Make a plot of the nodes
##'
##' Given the results of the splitting and testing algorithm in the form of a , make a node level data set for use in reporting results in terms of a binary tree graph. This does not print or plot the graph. You'll need to do that with the resulting object.
##' @param orig_res results data.table output from the \code{\link{findBlocks}} function.
##' @param fwer Means that alpha is fixed otherwise use the alpha that is calculated at each step of the splitting procedure
##' @return A ggraph object
make_graph <- function(res_graph) {
  require(ggraph)
  require(tidygraph)
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
