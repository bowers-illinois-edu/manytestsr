## Test how to report and plot. So far no actual tests here yet.

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(testthat)
  testthat::local_edition(3)
  library(here)
  library(data.table)
  library(dtplyr)
  library(dplyr)
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  setDTthreads(1)
  load_all() ## use  this during debugging
}

data(example_dat, package = "manytestsr")
example_dat$blockF <- factor(example_dat$blockF)
example_dat$trtF <- factor(example_dat$trtF)
example_dat <- droplevels(example_dat) %>% as.data.table()
## nrow(example_dat)

## Make the block-level dataset
example_bdat <- example_dat %>%
  group_by(blockF) %>%
  summarize(
    nb = n(),
    pb = mean(trt),
    hwt = (nb / nrow(example_dat)) * (pb * (1 - pb)),
    place = unique(place),
    year = unique(year),
    place_year_block = factor(unique(place_year_block))
  ) %>%
  as.data.table()
example_bdat

example_dat[, nb := .N, by = blockF]
table(example_dat$nb)

## dists1 <- example_dat[,dists_and_trans(Y1),by=blockF]
## dists2 <- example_dat[,fast_dists_and_trans_hybrid(Y1),by=blockF]
## expect_equal(dists1,dists2)
##
## uniq_ys <- example_dat %>% group_by(blockF) %>%
##  summarize(n_uniq_y1=length(unique(Y1))) %>%
##  arrange(n_uniq_y1) %>% as.data.table()
##
## dists <- cbind(d1=dists1,d2=dists2,example_dat[,.(Y1,nb)])
## all.equal(dists$d1.blockF,dists$d2.blockF)
## dists[,same_mean_rank_dist:=(d1.mean_rank_dist-d2.mean_rank_dist),by=d1.blockF]
##
## dists %>%
##  filter(d1.blockF %in% c("B092","B107","B108","B090")) %>%
##  select(Y1,nb,d1.rankY,d1.blockF,d1.mean_dist,d1.mean_rank_dist,d2.mean_rank_dist)
##
## summary(dists$same_mean_rank_dist,exclude=c())
## table(dists$same_mean_rank_dist==0)
## good_results <- dists[same_mean_rank_dist==0,.(d1.blockF,d1.rankY,Y1,nb,same_mean_rank_dist)]
## summary(good_results$nb)
## summary(good_results$d1.rankY)
##
## bad_results <- dists[same_mean_rank_dist!=0,.(d1.blockF,d1.rankY,Y1,nb,same_mean_rank_dist)]
## summary(bad_results$nb)
## summary(bad_results$d1.rankY)
##
##
## dists1 <- dists_and_trans(example_dat[nb>2,Y1])
## dists2 <- fast_dists_and_trans_hybrid(example_dat[nb>2,Y1])
## expect_equal(dists1,dists2)
##
##
## example_dat[,rank_Y1_v1:=Rfast::Rank(Y1),by=blockF]
## example_dat[,rank_Y1_v2:=avg_rank_arma(Y1),by=blockF]
##
## all.equal(example_dat$rank_Y1_v1, example_dat$rank_Y1_v2)
##
## table(example_bdat$nb)
##
## Varying block sizes, some variation in proportion in treatment within block
## Blocks created by place and year
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = fast_dists_and_trans_new)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = fast_dists_and_trans_hybrid)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = dists_and_trans)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF)
pIndepDist(dat = example_dat, fmla = Y2 ~ trtF | blockF)
## for comparison:
pOneway(dat = example_dat, fmla = Y1 ~ trtF | blockF)
pOneway(dat = example_dat, fmla = Y2 ~ trtF | blockF)

## Test the splitting functions
## splitSpecifiedFactor(example_bdat$blockF, example_bdat$place_year_block)
splitSpecifiedFactorMulti(example_bdat$blockF, example_bdat$place_year_block)
splitSpecified(example_bdat$blockF, example_bdat[, .(place, year, blockF)])

## find_blocks returns a data.table/data.frame with each row describing a block.

example_blocks_spec_fwer <- find_blocks(
  idat = example_dat, bdat = example_bdat, blockid = "blockF",
  pfn = pIndepDist, fmla = Y1 ~ trtF | blockF,
  splitfn = splitSpecifiedFactorMulti,
  alphafn = NULL,
  splitby = "place_year_block",
  blocksize = "hwt",
  copydts = TRUE, ncores = 1, parallel = "no", trace = FALSE
)

## Compare hits/discoveries from report_detections (which operates on blocks)
## and make_results_tree (which operates on a tree or graph object with nodes
## and leaves)

example_hits_spec_fwer <- report_detections(example_blocks_spec_fwer$bdat, fwer = TRUE, only_hits = FALSE, blockid = "blockF")
example_tree_spec_fwer <- make_results_tree(example_blocks_spec_fwer$bdat, block_id = "blockF")
make_results_ggraph(example_tree_spec_fwer$graph)
example_nodes_spec_fwer <- example_tree_spec_fwer$graph %>%
  activate(nodes) %>%
  as.data.frame()
## In this case we want to know the place and year and block info for nodes below the first overall test.
## That is, we learn (1) that we can reject the null of no effects overall
## (and that we have 5 splits from this first test)
## example_nodes_spec_fwer %>%
##  select(nodenum, name, parent_name, p, a, depth, out_degree, num_blocks, leaf_hit, group_hit, hit)
#### See just two nodes where we can reject the null: (node 1 --- the overall node, and node 6)
example_tree_spec_fwer$graph %>%
  activate(edges) %>%
  as.data.frame()

example_hits_spec_fwer[, .(hit, hit_grp, fin_grp, fin_nodenum, fin_parent, fin_parent_p, max_p)]

## Now, in this case we have one other test that rejects the null ---
## this time for the null of no effects among the units in the blocks under node 03c36456.
example_blocks_spec_fwer$bdat %>%
  filter(nodenum_prev == example_nodes_spec_fwer$name[6]) %>%
  select(nodenum_prev, place, year, blockF, pfinalb)

## We see that we don't have enough information to reject the null of no effects in the other places, but can reject the null of no effects for place A
example_blocks_spec_fwer$bdat %>%
  group_by(nodenum_prev) %>%
  reframe(unique(place), n())


## Compare to bottom-up fdr
### Notice that even with FDR we lose power. There are 2 (or 3) blocks with low
### unadjusted p-values but after adjustment we cannot reject

## example_blocks_fdr <- adjust_block_tests(
##  idat = example_dat, bdat = example_bdat, blockid = "blockF",
##  pfn = pIndepDist,
##  p_adj_method = "fdr",
##  fmla = Y1 ~ trtF,
##  copydts = TRUE, ncores = 1, parallel = "no"
## )
## example_blocks_fdr %>%
##  select(p, max_p) %>%
##  arrange(p) %>%
##  head()

example_blocks_spec_fdr <- find_blocks(
  idat = example_dat, bdat = example_bdat, blockid = "blockF",
  pfn = pIndepDist, fmla = Y1 ~ trtF | blockF,
  splitfn = splitSpecifiedFactorMulti,
  alphafn = alpha_saffron,
  splitby = "place_year_block",
  blocksize = "nb", thealpha = .05, thew0 = .05 - .00000001,
  copydts = TRUE, ncores = 1, parallel = "no", trace = FALSE
)$bdat
example_hits_spec_fdr <- report_detections(example_blocks_spec_fdr, fwer = FALSE, alpha = .1, only_hits = FALSE, blockid = "blockF")

example_tree_spec_fdr <- make_results_tree(example_blocks_spec_fdr, block_id = "blockF")
make_results_ggraph(example_tree_spec_fdr$graph)
example_nodes_spec_fdr <- example_tree_spec_fdr$graph %>%
  activate(nodes) %>%
  as.data.frame()
## In this case we want to know the place and year and block info for nodes below the first overall test.
## That is, we learn (1) that we can reject the null of no effects overall
## (and that we have 5 splits from this first test)
example_nodes_spec_fdr %>%
  select(node_number, name, parent_name, p, a, depth, num_leaves)
## See just two nodes where we can reject the null: (node 1 --- the overall node, and node 6)
example_tree_spec_fdr$graph %>%
  activate(edges) %>%
  as.data.frame()

node_info_leaves <- example_hits_spec_fdr[, .(
  nodesize = sum(unique(nodesize)),
  places = paste(unique(place), collapse = ",")
), by = .("nodenum" = nodenum_prev)]
node_info_parent <- example_hits_spec_fdr[, .(
  nodesize = sum(unique(nodesize)),
  places = paste(unique(place), collapse = ",")
), by = .("nodenum" = nodenum_current)]
node_info <- rbind(node_info_parent, node_info_leaves)
stopifnot(length(unique(node_info$nodenum)) == nrow(node_info))

# example_nodes_spec_fdr2 <- merge(example_nodes_spec_fdr, node_info, by = "node_number", all = TRUE)
# example_nodes_spec_fdr2$nodesize[example_nodes_spec_fdr2$nodenum == "1"] <- nrow(example_dat)
# example_nodes_spec_fdr2 %>%
#  select(nodenum, name, parent_name, p, a, depth, out_degree, nodesize, places) %>%
#  arrange(depth)
#
#
# with(
#  example_nodes_spec_fdr2,
#  alpha_saffron(pval = p, batch = depth, nodesize = nodesize, thealpha = .1, thew0 = .1 - .001)
# )
#

# ggraph(example_tree_spec_fdr$graph) +
#  geom_edge_diagonal() +
#  geom_node_label(aes(label = round(p, 3), colour = hit),
#    repel = FALSE, show.legend = FALSE, label.r = unit(0.5, "lines"),
#    label.padding = unit(.01, "lines"), label.size = 0
#  )
###


### TODO: Test the error and rejection and power calculations
