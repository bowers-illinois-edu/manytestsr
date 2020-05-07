
break


# res3g <- ToDiagrammeRGraph(res3tree2)
# orig_nodes_df <- res3g %>% get_node_df()
# new_nodes_df <- merge(orig_nodes_df, nodes_df, by = "label", all.x = TRUE, sort = FALSE)
# new_nodes_df <- new_nodes_df[, c(names(orig_nodes_df), "p", "bF","depth")]
# stopifnot(all.equal(new_nodes_df[, 1:3], orig_nodes_df))
# res3g$nodes_df <- new_nodes_df
# res3g$nodes_df$label <- paste(res3g$nodes_df$label, "\n Blocks:", res3g$nodes_df$bF, "\n p=", round(res3g$nodes_df$p, 4), sep = "")
## res3g <- res3g %>% set_node_attrs(node_attr=color,values="green",nodes=which(res3g$nodes_df$p < .05))
# res3g$nodes_df$type <- ifelse(res3g$nodes_df$p < .05, "NS", "S")
# res3g$nodes_df$fillcolor <- ifelse(res3g$nodes_df$p < .05, gray(.5), gray(.8))
# res3g$nodes_df$style <- rep("filled", count_nodes(res3g))
# res3g$nodes_df$fontcolor <- "black"
# res3g <- res3g %>%
#   select_nodes_by_degree(expressions = "deg==1") %>%
#   select_nodes(conditions = p <= .05, set_op = "intersect") %>%
#   set_node_attrs_ws(node_attr = fillcolor, value = "orange")
## If both leaves are p>alpha BUT their parent is p<alpha,then report BOTH of them as a set.
# res3g <- res3g %>% join_node_attrs(df = data.frame(get_degree_in(.),
#                                   get_degree_out(.)))
# render_graph(res3g)
#
# library(igraph)
# library(ggraph)
# library(tidygraph)
# library(ggtree)
# res3i <- to_igraph(res3g)

res3gg <- as_tbl_graph(res3i)

res3gg %>% mutate(root = node_is_root(), leaf = node_is_leaf())


res3gg %>%
  activate(nodes) %>%
  mutate(P = map_bfs_chr(node_is_root(),
    .f = function(node, path, ...) {
      .N()$Employee[tail(path$node[.N()$Grade[node] - .N()$Grade[path$node] >= 2], 1)[1]]
    }
  ))

g <- ggraph(res3gg, layout = "tree") +
  geom_edge_link() +
  # geom_node_point() +
  geom_node_label(aes(label = label))
print(g)

gt <- ggtree(res3i, layout = "slanted")
print(gt)


# sample graph
g <- tidygraph::create_tree(8, 2) %>%
  activate(nodes) %>%
  mutate(
    Employee = LETTERS[1:8],
    Grade = c(0, 1, 1, 1, 2, 2, 2, 2),
    label = paste("Emp", Employee, "Grade", Grade)
  )
library(visNetwork)
visIgraph(g, layout = "layout_as_tree", flip.y = F, idToLabel = F)

g %>%
  activate(nodes) %>%
  mutate(P = map_bfs_chr(node_is_root(), .f = function(node, path, ...) {
    .N()$Employee[tail(path$node[.N()$Grade[node] - .N()$Grade[path$node] >= 2], 1)[1]]
  })) %>%
  as_tibble()

plot(res3i,
  layout = layout_as_tree,
  vertex.shape = "rectangle",
  vertex.size = (strwidth(V(res3i)$label) + strwidth("oo")) * 100,
  vertex.size2 = strheight(V(res3i)$label) * 2 * 100
)


break

bdatTrueATEs <- idat3[, lapply(.SD, mean), by = bF, .SDcols = grep("^tau", names(idat3), value = TRUE)]

test_that("alphafns work across splitters", {
  resnms <- apply(alpha_and_splits, 1, function(x) {
    paste(x, collapse = "_", sep = "")
  })
  ## First with no effects at all. So, basically no splitting should happen.
  bigres_null <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, fmla = Ynull ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(bigres_null) <- resnms

  bigres_homog <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, fmla = Yhomog ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(bigres_homog) <- resnms

  bigres_normb <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, fmla = Ynormb ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(bigres_normb) <- resnms
})
#####
# res3$maxpC <- round(res3$maxp, 3)
# res3tree <- as.Node(res3[, .(bF, maxp, biggrp)], pathDelimiter = ".", pathName = "biggrp")
res3tree2 <- as.Node(res3long, pathDelimiter = ".", pathName = "biggrp")
print(res3tree2, "a", "p", "depth", "bFC")
res3tree2$Get("p")

res3i <- as.igraph(res3tree2, vertexAttributes = c("p", "depth", "bfC", "a"), directed = TRUE)
res3g <- as_tbl_graph(res3i)
blahV <- res3g %>%
  activate(nodes) %>%
  as_tibble()
blahE <- res3g %>%
  activate(edges) %>%
  as_tibble()
res3g <- res3g %>%
  activate(nodes) %>%
  mutate(label = paste(name, "\n Blocks:", bfC, "\n p=", round(p, 4), sep = ""))


g <- ggraph(res3g, layout = "tree") +
  geom_edge_link() +
  # geom_node_point() +
  geom_node_label(aes(label = label))
print(g)


##
## res3nodes_lst <- list()
## for (i in 1:max(res3long$depth)) {
##   nodenm <- paste("nodenum", i, sep = "")
##   res3nodes_lst[[i]] <- res3long[!is.na(get(nodenm)) & depth == i,
##     .(
##       p = unique(p),
##       bF = paste(as.character(unlist(sort(bF))), collapse = ","),
##       depth = unique(depth)
##     ),
##     by = get(nodenm)
##   ]
## }
## res3nodes_df <- rbindlist(res3nodes_lst)
## setnames(res3nodes_df, "get", "labelN")
## res3nodes_df$label <- as.character(res3nodes_df$labelN)
## res3nodes_df$name <- res3nodes_df$label
##
