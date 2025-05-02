## Test how to report and plot. So far no actual tests here yet.
context("Reporting and Plotting")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
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
example_dat <- droplevels(example_dat)
## nrow(example_dat)

## Make the block-level dataset
example_bdat <- example_dat %>%
  group_by(blockF) %>%
  summarize(
    nb = n(),
    pb = mean(trt),
    hwt = (nb / nrow(example_dat)) * (pb * (1 - pb)),
    site = unique(site),
    year = unique(year),
    site_year_block = unique(site_year_block)
  ) %>%
  as.data.table()
example_bdat

## Varying block sizes, some variation in proportion in treatment within block
## Blocks created by site and year
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = fast_dists_and_trans_by_unit_arma2, adaptive_dist_function = FALSE)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = dists_and_trans, adaptive_dist_function = FALSE)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF, distfn = fast_dists_and_trans_by_unit_arma_parR(threads = 4), adaptive_dist_function = FALSE)
pIndepDist(dat = example_dat, fmla = Y1 ~ trtF | blockF)
pIndepDist(dat = example_dat, fmla = Y2 ~ trtF | blockF)
## for comparison:
pOneway(dat = example_dat, fmla = Y1 ~ trtF | blockF)
pOneway(dat = example_dat, fmla = Y2 ~ trtF | blockF)

## Test the splitting functions
## splitSpecifiedFactor(example_bdat$blockF, example_bdat$site_year_block)
splitSpecifiedFactorMulti(example_bdat$blockF, example_bdat$site_year_block)
splitSpecified(example_bdat$blockF, example_bdat[, .(site, year, blockF)])

## find_blocks returns a data.table/data.frame with each row describing a block.

example_blocks_spec_fwer <- find_blocks(
  idat = example_dat, bdat = example_bdat, blockid = "blockF",
  pfn = pIndepDist, fmla = Y1 ~ trtF | blockF,
  splitfn = splitSpecifiedFactorMulti,
  alphafn = NULL,
  splitby = "site_year_block",
  blocksize = "hwt",
  copydts = TRUE, ncores = 1, parallel = "no", trace = FALSE
)

## Compare hits/discoveries from report_detections (which operates on blocks)
## and make_results_tree (which operates on a tree or graph object with nodes
## and leaves)

example_hits_spec_fwer <- report_detections(example_blocks_spec_fwer, fwer = TRUE, only_hits = FALSE, blockid = "blockF")
example_tree_spec_fwer <- make_results_tree(example_blocks_spec_fwer, blockid = "blockF")
make_results_ggraph(example_tree_spec_fwer)
example_nodes_spec_fwer <- example_tree_spec_fwer %>%
  activate(nodes) %>%
  as.data.frame()
## In this case we want to know the site and year and block info for nodes below the first overall test.
## That is, we learn (1) that we can reject the null of no effects overall
## (and that we have 5 splits from this first test)
example_nodes_spec_fwer %>%
  select(nodenum, name, parent_name, p, a, depth, out_degree, num_blocks, leaf_hit, group_hit, hit)
## See just two nodes where we can reject the null: (node 1 --- the overall node, and node 6)
example_tree_spec_fwer %>%
  activate(edges) %>%
  as.data.frame()

example_hits_spec_fwer[, .(hit, hit_grp, fin_grp, fin_nodenum, fin_parent, fin_parent_p, max_p)]

## Now, in this case we have one other test that rejects the null ---
## this time for the null of no effects among the units in the blocks under node 03c36456.
example_blocks_spec_fwer %>%
  filter(nodenum_prev == example_nodes_spec_fwer$nodenum[6]) %>%
  select(nodenum_prev, site, year, blockF, pfinalb)

## We see that we don't have enough information to reject the null of no effects in the other sites, but can reject the null of no effects for site A
example_blocks_spec_fwer %>%
  group_by(nodenum_prev) %>%
  reframe(unique(site), n())


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
  splitby = "site_year_block",
  blocksize = "nb", thealpha = .05, thew0 = .05 - .00000001,
  copydts = TRUE, ncores = 1, parallel = "no", trace = FALSE
)
example_hits_spec_fdr <- report_detections(example_blocks_spec_fdr, fwer = FALSE, alpha = .1, only_hits = FALSE, blockid = "blockF")

example_tree_spec_fdr <- make_results_tree(example_blocks_spec_fdr, blockid = "blockF")
make_results_ggraph(example_tree_spec_fdr)
example_nodes_spec_fdr <- example_tree_spec_fdr %>%
  activate(nodes) %>%
  as.data.frame()
## In this case we want to know the site and year and block info for nodes below the first overall test.
## That is, we learn (1) that we can reject the null of no effects overall
## (and that we have 5 splits from this first test)
example_nodes_spec_fdr %>%
  select(nodenum, name, parent_name, p, a, depth, out_degree)
## See just two nodes where we can reject the null: (node 1 --- the overall node, and node 6)
example_tree_spec_fdr %>%
  activate(edges) %>%
  as.data.frame()

node_info_leaves <- example_hits_spec_fdr[, .(
  nodesize = sum(unique(nodesize)),
  sites = paste(unique(site), collapse = ",")
), by = .("nodenum" = nodenum_prev)]
node_info_parent <- example_hits_spec_fdr[, .(
  nodesize = sum(unique(nodesize)),
  sites = paste(unique(site), collapse = ",")
), by = .("nodenum" = nodenum_current)]
node_info <- rbind(node_info_parent, node_info_leaves)
stopifnot(length(unique(node_info$nodenum)) == nrow(node_info))

example_nodes_spec_fdr2 <- merge(example_nodes_spec_fdr, node_info, by = "nodenum", all = TRUE)
example_nodes_spec_fdr2$nodesize[example_nodes_spec_fdr2$nodenum == "1"] <- nrow(example_dat)
example_nodes_spec_fdr2 %>%
  select(nodenum, name, parent_name, p, a, depth, out_degree, nodesize, sites) %>%
  arrange(depth)


with(
  example_nodes_spec_fdr2,
  alpha_saffron(pval = p, batch = depth, nodesize = nodesize, thealpha = .1, thew0 = .1 - .001)
)


ggraph(example_tree_spec_fdr) +
  geom_edge_diagonal() +
  geom_node_label(aes(label = round(p, 3), colour = hit),
    repel = FALSE, show.legend = FALSE, label.r = unit(0.5, "lines"),
    label.padding = unit(.01, "lines"), label.size = 0
  )


## Setup the fake data to explore other forms of reporting
## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]

## Not sure what kind of test to write here. Leaving this as is for now. Hoping it doesn't break package building
bdat4[, x1 := 1:nrow(bdat4)]
bdat4[, x2 := c(rep(1, nrow(bdat4) - 5), rep(10, 5))]
bdat4[, x3 := c(rep(1, nrow(bdat4) - 5), 5, rep(10, 4))]
bdat4[, g5 := splitCluster(bid = as.character(bF), x = v4)]
bdat4[, g_x1 := splitCluster(bid = as.character(bF), x = x1)]
bdat4[, g_x2 := splitCluster(bid = as.character(bF), x = x2)]
bdat4[, g_x3 := splitCluster(bid = as.character(bF), x = x3)]
## Setting up  a test of pre-specified splits
## First, where each level is identified by its own column.
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
## This next for the splitSpecifiedFactor
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]
