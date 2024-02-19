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
  # conflicts_prefer(dplyr::filter)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  setDTthreads(1)
  load_all() ## use  this during debugging
}
library(here)

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

with(bdat4, table(lv1, lv3))
with(bdat4, table(lv1, lv2))
table(bdat4$lv1)

## Test pre-specified splitters
bdat4[, gf5 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf5))
bdat4[lv1 != "l1_1", gf6 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf6, exclude = c()))
bdat4[lv1 != "l1_1" & lv2 != "l2_2", gf7 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf7, exclude = c()))

bdat4[, gf8 := splitSpecifiedFactor(bF, x = lvs)]
with(bdat4, table(lvs, gf8))
bdat4[lv1 != "l1_1", gf9 := splitSpecifiedFactor(bF, x = lvs)]
with(bdat4, table(lvs, gf9, exclude = c()))
bdat4[lv1 != "l1_1" & lv2 != "l2_2", gf10 := splitSpecifiedFactor(bF, x = lvs)]
with(bdat4, table(lvs, gf10, exclude = c()))

## A splitting variable with little variation.
set.seed(12345)
bdat4[, twosplits := rbinom(.N, 1, .5)]
bdat4[, twosplitsF := factor(twosplits)] ## should only have 2 splits
bdat4[, lvs2 := interaction(lv1, lv2)] ## should only have 4 spits
bdat4[, constv := rep(10, .N)]


## Not testing "splitSpecified" because splitSpecifiedFactor does a better job.
## depreciate splitSpecified
split_test_parms <- data.table(expand.grid(
  sfn = c(
    "splitLOO", "splitEqualApprox",
    "splitCluster", "splitSpecifiedFactor"
  ),
  splitby = c("twosplits", "twosplitsF", "lvs2", "constv"),
  stopsplitting = c(TRUE, FALSE), stringsAsFactors = FALSE
))

split_test_parms <- split_test_parms[sfn != "splitSpecifiedFactor" |
  (sfn == "splitSpecifiedFactor" & !(splitby %in% c("twosplits", "constv"))), ]

test_splitters_fn <- function(sfn, splitby, stopsplitting) {
  theres <- findBlocks(
    idat = idat3, bdat = bdat4, blockid = "bF",
    pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla = Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE,
    splitfn = get(sfn), splitby = splitby, stop_splitby_constant = stopsplitting
  )
  return(theres)
}

res <- mapply(
  FUN = function(sfn = sfn, sby = sby, stopsplitting = stopsplitting) {
    message(paste(sfn, sby, stopsplitting, collapse = ","))
    ## Some errors are expected here.
    obj <- try(test_splitters_fn(
      sfn = sfn,
      splitby = sby,
      stopsplitting = stopsplitting
    ), silent = TRUE)
    obj_det <- if (class(obj)[1] == "try-error") {
      return(NA)
    } else {
      return(report_detections(obj))
    }
  },
  sfn = split_test_parms$sfn,
  sby = split_test_parms$splitby,
  stopsplitting = split_test_parms$stopsplitting,
  SIMPLIFY = FALSE
)

names(res) <- apply(split_test_parms, 1, function(x) {
  paste(x, collapse = "_", sep = "")
})

tree1 <- make_tree(res[["splitSpecifiedFactor_lvs2_ TRUE"]], blockid = "bF")
det1 <- report_detections(res[["splitSpecifiedFactor_lvs2_FALSE"]], blockid = "bF", only_hits = TRUE, fwer = TRUE, alpha = .05)

tree1tib <- tree1 %>%
  activate(nodes) %>%
  as_tibble()

load(here("tests", "dpp_dat.rda"))
dpp_dat <- droplevels(dpp_dat)

dpp_bdat <- dpp_dat %>%
  group_by(blockF) %>%
  summarize(
    nb = n(),
    pb = mean(trt),
    hwt = (nb / nrow(dpp_dat)) * (pb * (1 - pb)),
    site = unique(STSITE),
    cohort = unique(STCOHORT),
    site_cohort_block = unique(site_cohort_block)
  ) %>%
  as.data.table()
head(dpp_bdat)

pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
pIndepDist(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)

dpp_blocks_1 <- findBlocks(
  idat = dpp_dat, bdat = dpp_bdat, blockid = "blockF",
  pfn = pWilcox, fmla = R01TMCRET ~ trtF | blockF,
  splitfn = splitSpecifiedFactor,
  alphafn = NULL,
  splitby = "site_cohort_block",
  blocksize = "hwt",
  copydts = TRUE, ncores = 8, parallel = "multicore", trace = FALSE
)


detections_dpp <- report_detections(dpp_blocks_1, fwer = TRUE, only_hits = FALSE, blockid = "blockF")

tmp <- detections_dpp %>%
  select(blockF, site_cohort_block, pfinalb, fin_nodenum, fin_parent, fin_parent_p, max_p, hit, fin_grp, hit_grp) %>%
  arrange(fin_nodenum)

dpp_blocks_1_tree <- make_tree(dpp_blocks_1, blockid = "blockF")


make_graph(make_tree(dpp_blocks_1, blockid = "blockF"))

tmp2 <- dpp_blocks_1_tree %>%
  activate(nodes) %>%
  as_tibble()
tmp3 <- left_join(tmp2, dpp_bdat[, c("blockF", "site", "cohort", "site_cohort_block")], by = join_by(bF == blockF))

dpp_blocks_1_tree <- dpp_blocks_1_tree %>%
  activate(nodes) %>%
  mutate(bFnew = stri_replace_all(bF, "", regex = "Block0+"))

ggraph(dpp_blocks_1_tree, layout = "tree") +
  geom_edge_diagonal() +
  geom_node_label(aes(label = round(p, 3), colour = hit),
    repel = FALSE, show.legend = FALSE, label.r = unit(0.5, "lines"),
    label.padding = unit(.01, "lines"), label.size = 0
  )

ggraph(dpp_blocks_1_tree, layout = "dendrogram") +
  geom_edge_diagonal() +
  geom_node_point(aes(filter = leaf)) +
  coord_fixed()
