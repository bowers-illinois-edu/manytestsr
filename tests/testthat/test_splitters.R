# Test and develop functions to split the data
testthat::context("Performance of Splitting Functions")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(dtplyr)
  library(dplyr)
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  setDTthreads(1)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  load_all() ## use  this during debugging
}
library(here)
library(dtplyr)
library(dplyr)

## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]

test_that("Leave One Outsplitting actually Leaves One out.", {
  ## The test of leave one out splitting is easiest.
  test_Split_LOO <- splitLOO(bdat4$bF, bdat4$hwt)
  expect_equal(sum(test_Split_LOO == "1"), length(test_Split_LOO) - 1)
})

## Not sure what kind of test to write here. Leaving this as is for now.
##  Hoping it doesn't break package building
bdat4[, x1 := 1:nrow(bdat4)]
bdat4[, x2 := c(rep(1, nrow(bdat4) - 5), rep(10, 5))]
bdat4[, x3 := c(rep(1, nrow(bdat4) - 5), 5, rep(10, 4))]
bdat4[, g5 := splitCluster(bid = as.character(bF), x = v4)]
bdat4[, g_x1 := splitCluster(bid = as.character(bF), x = x1)]
bdat4[, g_x2 := splitCluster(bid = as.character(bF), x = x2)]
bdat4[, g_x3 := splitCluster(bid = as.character(bF), x = x3)]
with(bdat4, table(g_x1, x1, exclude = c()))
with(bdat4, table(g_x2, x2, exclude = c()))
with(bdat4, table(g_x3, x3, exclude = c()))


test_that("Clustering based splitting makes the splits more homogeneous and distant in mean", {
  bdat4[, g1 := splitCluster(bid = as.character(bF), x = hwt)]
  bdat4[, g2 := splitEqualApprox(bid = as.character(bF), x = hwt)]
  # Maybe we want something from splitCLuster like the difference between the  means is bigger than any random split?
  g1summary <- bdat4[, .(mn = mean(hwt), sd = sd(hwt)), by = g1]
  g2summary <- bdat4[, .(mn = mean(hwt), sd = sd(hwt)), by = g2]
  expect_lt(max(g1summary$sd), min(g2summary$sd))
  expect_gt(abs(diff(g1summary$mn)), abs(diff(g2summary$mn)))
})


## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]

with(bdat4, table(lv1, lv3))
with(bdat4, table(lv1, lv2))
table(bdat4$lv1)
table(bdat4$lvs)

## Test pre-specified splitters
## deprecate splitSpecified
bdat4[, gf5 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf5))
bdat4[lv1 != "l1_1", gf6 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf6, exclude = c()))
bdat4[lv1 != "l1_1" & lv2 != "l2_2", gf7 := splitSpecified(bF, x = data.table(lv1, lv2, bF))]
with(bdat4, table(lvs, gf7, exclude = c()))

## The first split should be from the left-most element of the label
bdat4[, gf8 := splitSpecifiedFactorMulti(bF, x = lvs)]
with(bdat4, table(lvs, gf8))
bdat4[lv1 != "l1_1", gf9 := splitSpecifiedFactorMulti(bF, x = lvs)]
with(bdat4, table(lvs, gf9, exclude = c()))
with(droplevels(bdat4[lv1 != "l1_1", ]), table(lvs, gf9))
bdat4[lv1 != "l1_1" & lv2 != "l2_2", gf10 := splitSpecifiedFactorMulti(bF, x = lvs)]
with(bdat4, table(lvs, gf10, exclude = c()))
## This should have five splits
with(droplevels(bdat4[lv1 != "l1_1" & lv2 != "l2_2", ]), table(lvs, gf10, exclude = c()))

data(example_dat, package = "manytestsr")
example_dat$blockF <- factor(example_dat$blockF)

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

example_bdat[, g1 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(example_bdat, table(g1, site, exclude = c()))
# with(example_bdat, table(g1, site_year_block, exclude = c()))
example_bdat[g1 == 0, g2 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(example_bdat, table(g2, site, exclude = c()))
example_bdat[g2 == 0, g3 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(example_bdat, table(g3, site, exclude = c()))
example_bdat[g3 == 0, g4 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(example_bdat, table(g4, site, exclude = c()))
example_bdat[g4 == 0, g5 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(droplevels(example_bdat[g4 == 0, ]), table(g5, site_year_block, exclude = c()))
example_bdat[g5 == 0, g6 := splitSpecifiedFactor(bid = blockF, x = site_year_block)]
with(droplevels(example_bdat[g5 == 0, ]), table(g6, site_year_block, exclude = c()))

#### When splitby has no variation within biggrp, further splitting doesn't do
### anything and the algorithm may just continue to produce the same p-values
### until maxtest is reached depending on the splitting function. So, for
### example, in splitSpecifiedFactor and splitCluster we want the algorithm to
### stop or switch to another splitter (TODO maybe) once clusters / branches
### have no variation on splitby  --- any further splitting would be essentially
### random.

## In splitSpecifiedFactor we want it to stop when splitby is constant within group
## In splitCluster we want it to stop when splitby is constant within group or move
## to random splits or splitLOO or splitEqualApprox
## In splitEqualApprox we want splitby to be something like N and so we want it to
## continue even if the groups have the same N.

## In splitLOO could stop when splitby is constant within the bigger group or
## start choosing blocks to focus on at random.
## for the first split (say, with only two groups or even with only one group),
## splitLOO will choose one block at random. So, splitLOO with a splitting
## vector/criteria/variable that does not vary is the same as random splits.
## So, probably ok for error rate but not great for power.

## splitting variables with little variation.
set.seed(12345)
bdat4[, twosplits := rbinom(.N, 1, .5)]
bdat4[, twosplitsF := factor(twosplits)] ## should only have 2 splits
bdat4[, lvs2 := interaction(lv1, lv2)] ## should only have 4 spits
bdat4[, constv := rep(10, .N)] ## a constant

test_that("splitCluster follows the values of discrete splitby variables", {
  theres1 <- findBlocks(
    idat = idat3, bdat = bdat4, blockid = "bF", pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla = Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE,
    splitfn = splitCluster, splitby = "twosplits", stop_splitby_constant = TRUE
  )
  theres1_det <- report_detections(theres1, blockid = "bF")
  expect_equal(uniqueN(theres1$biggrp), uniqueN(bdat4$twosplits))
  ## This next because I know that we should reject in both groups.
  expect_equal(uniqueN(theres1_det$fin_grp), uniqueN(bdat4$twosplits))
  thetab <- table(theres1$twosplitsF, theres1$biggrp)
  expect_equal(c(thetab[1, 2], thetab[2, 1]), c(0, 0))
})


## This next doesn't have explicit expectations and so is a skipped or empty test upon R check. Comment out for now.
test_that("splitCluster stops appropriately (i.e. doesn't just keep randomly splitting) with continuous splitting criteria.", {
  theres1 <- findBlocks(
    idat = idat3, bdat = bdat4, blockid = "bF", pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla = Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE,
    splitfn = splitCluster, splitby = "hwt", stop_splitby_constant = TRUE
  )
  theres1_det <- report_detections(theres1, blockid = "bF")
  theres2 <- findBlocks(
    idat = idat3, bdat = bdat4, blockid = "bF", pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla = Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE,
    splitfn = splitCluster, splitby = "hwt", stop_splitby_constant = FALSE
  )
  theres2_det <- report_detections(theres1, blockid = "bF")
  ##   ## Not sure what else to expect here
  expect_equal(theres1_det, theres2_det)
})



## Not testing "splitSpecified" because splitSpecifiedFactor and splitSpecifiedFactorMulti do better jobs.
## depreciate splitSpecified
split_test_params <- data.table(expand.grid(
  sfn = c(
    "splitLOO", "splitEqualApprox",
    "splitCluster", "splitSpecifiedFactor",
    "splitSpecifiedFactorMulti"
  ),
  splitby = c("twosplits", "twosplitsF", "lvs2", "constv"),
  stopsplitting = c(TRUE, FALSE), stringsAsFactors = FALSE
))

split_test_params <- split_test_params[sfn != "splitSpecifiedFactorMulti" |
  (sfn == "splitSpecifiedFactorMulti" & !(splitby %in% c("twosplits", "constv"))), ]

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

res1 <- test_splitters_fn(
  sfn = split_test_params$sfn[1],
  splitby = split_test_params$splitby[1],
  stopsplitting = split_test_params$stopsplitting[1]
)
res1_det <- report_detections(res1, blockid = "bF")

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
  sfn = split_test_params$sfn,
  sby = split_test_params$splitby,
  stopsplitting = split_test_params$stopsplitting,
  SIMPLIFY = FALSE
)

names(res) <- apply(split_test_params, 1, function(x) {
  gsub(" ", "", paste0(x, collapse = "_")) ## DK why adding space.
})

test_that("Splitters work as expected given splitby variables
  that have one, two or four values and are factor or numeric", {
  ### Expectations ###
  ## > split_test_params
  ##                      sfn    splitby stopsplitting
  ## splitLOO will tend to first choose one of twosplits==1 at random and then
  ## again until there are no more twosplits==1, then it will stop. If it would
  ## reject with the entire group of twosplits==0, it will still stop but then
  ## report that as a rejection.
  ##  1:             splitLOO  twosplits          TRUE
  resobj <- which(names(res) == "splitLOO_twosplits_TRUE")
  with(res[[resobj]], table(biggrp, twosplits, exclude = c()))
  with(res[[resobj]], table(fin_grp, twosplits, exclude = c()))
  ## All of the group with the smallest value on twosplits will be put in one split
  expect_equal(sum(res[[resobj]]$twosplits == 0), max(table(res[[resobj]]$fin_grp)))
  ##  4:             splitLOO twosplitsF          TRUE
  ### This should be the same as above. Factor is the same as numeric here.
  resobj <- which(names(res) == "splitLOO_twosplits_TRUE")
  with(res[[resobj]], table(biggrp, twosplits, exclude = c()))
  with(res[[resobj]], table(fin_grp, twosplits, exclude = c()))
  expect_equal(sum(res[[resobj]]$twosplitsF == 0), max(table(res[[resobj]]$fin_grp)))

  ## 15:             splitLOO  twosplits         FALSE
  ### This should be like random splits: Yes. See below how blah2 and blah2a differ.
  ## So, splitLOO really should be used with a splitby that is more continuous and random splits will
  ## happen when the systematic variation in splitby is used up.
  resobj <- which(names(res) == "splitLOO_twosplits_FALSE")
  with(res[[resobj]], table(biggrp, twosplits, exclude = c()))
  with(res[[resobj]], table(fin_grp, twosplits, exclude = c()))
  with(res[[resobj]], table(hit, twosplits, exclude = c()))

  blah2 <- test_splitters_fn(sfn = "splitLOO", splitby = "twosplits", stopsplitting = FALSE)
  blah2_det <- report_detections(blah2, blockid = "bF")
  with(blah2_det, table(fin_grp, twosplits, exclude = c()))

  blah2a <- test_splitters_fn(sfn = "splitLOO", splitby = "twosplits", stopsplitting = FALSE)
  blah2a_det <- report_detections(blah2a, blockid = "bF")
  with(blah2a_det, table(fin_grp, twosplits, exclude = c()))
  blah2_tree <- make_graph(make_tree(blah2))


  ## 18:             splitLOO twosplitsF         FALSE

  resobj <- which(names(res) == "splitLOO_twosplitsF_FALSE")
  with(res[[resobj]], table(biggrp, twosplitsF, exclude = c()))
  with(res[[resobj]], table(fin_grp, twosplitsF, exclude = c()))
  with(res[[resobj]], table(hit, twosplitsF, exclude = c()))

  blah2F <- test_splitters_fn(sfn = "splitLOO", splitby = "twosplitsF", stopsplitting = FALSE)
  blah2F_det <- report_detections(blah2F, blockid = "bF")
  with(blah2F_det, table(fin_grp, twosplits, exclude = c()))

  blah2Fa <- test_splitters_fn(sfn = "splitLOO", splitby = "twosplitsF", stopsplitting = FALSE)
  blah2Fa_det <- report_detections(blah2Fa, blockid = "bF")
  with(blah2Fa_det, table(fin_grp, twosplits, exclude = c()))

  ##  8:             splitLOO       lvs2          TRUE
  resobj <- which(names(res) == "splitLOO_lvs2_TRUE")
  ### Basically don't use splitLOO with overly discrete splitby unless we are happy with random splits
  with(res[[resobj]], table(biggrp, lvs2, exclude = c()))
  with(res[[resobj]], table(fin_grp, lvs2, exclude = c()))

  ## 12:             splitLOO     constv          TRUE
  ### This one failed because you can't have a splitby that is constant if stop_splitby_constant=TRUE
  ## blah2G <- test_splitters_fn(sfn = "splitLOO", splitby = "constv", stopsplitting = TRUE)
  expect_equal(res[["splitLOO_constv_TRUE"]], NA)

  ## 22:             splitLOO       lvs2         FALSE
  ## Semi-random splits after exhausting variation in lvs2
  with(res[["splitLOO_lvs2_FALSE"]], table(fin_grp, lvs2, exclude = c()))

  ## 26:             splitLOO     constv         FALSE
  ## This is random splits again.
  with(res[["splitLOO_constv_FALSE"]], table(fin_grp, constv, exclude = c()))
  ###########

  ############################

  ### splitEqualApprox tries to make two groups such that the sum(x_{group
  ### 1})=sum(x_{group_2}). This means that in a case with only two values, it
  ### tries to allocate blocks to splits such that the numbers of blocks with x=0
  ### and x=1 in group 1 are the same as they are in group 2
  ### Unclear what patterns to expect here
  ##  2:     splitEqualApprox  twosplits          TRUE
  with(res[["splitEqualApprox_twosplits_TRUE"]], table(fin_grp, twosplits, exclude = c()))
  blah2H <- test_splitters_fn(
    sfn = "splitEqualApprox", splitby = "twosplits",
    stopsplitting = TRUE
  )
  blah2H_det <- report_detections(blah2H, blockid = "bF")
  with(blah2H_det, table(fin_grp, twosplits, exclude = c()))

  ##  5:     splitEqualApprox twosplitsF          TRUE
  with(res[["splitEqualApprox_twosplitsF_TRUE"]], table(fin_grp, twosplitsF, exclude = c()))

  ##  9:     splitEqualApprox       lvs2          TRUE
  ### splitEqualApprox will rank factors and sum them too. So, again, not really that great for factors.
  with(res[["splitEqualApprox_lvs2_TRUE"]], table(fin_grp, lvs2, exclude = c()))

  ## 13:     splitEqualApprox     constv          TRUE
  ## Not allowed by findBlocks.
  expect_equal(res[["splitEqualApprox_constv_TRUE"]], NA)

  ## 16:     splitEqualApprox  twosplits         FALSE
  ## Random splits basically  (after first split then random)
  with(res[["splitEqualApprox_twosplits_FALSE"]], table(fin_grp, twosplits, exclude = c()))

  ## 19:     splitEqualApprox twosplitsF         FALSE
  ### Random splits again (after first split then random)
  with(res[["splitEqualApprox_twosplitsF_FALSE"]], table(fin_grp, twosplits, exclude = c()))

  ## 23:     splitEqualApprox       lvs2         FALSE
  ## random splits after initial splits by lvs2 (ranking levels,
  ## taking every other level into group 1 the rest into group 2)
  with(res[["splitEqualApprox_lvs2_FALSE"]], table(fin_grp, lvs2, exclude = c()))

  ## 27:     splitEqualApprox     constv         FALSE
  ### random splits
  with(res[["splitEqualApprox_constv_FALSE"]], table(fin_grp, constv, exclude = c()))

  ##########################
  ##  3:         splitCluster  twosplits          TRUE
  tabres3 <- with(res[["splitCluster_twosplits_TRUE"]], table(fin_grp, twosplits, exclude = c()))
  expect_equal(tabres3[1, 2], 0)
  expect_equal(tabres3[2, 1], 0)
  expect_equal(dim(tabres3), c(2, 2))

  ##  6:         splitCluster twosplitsF          TRUE
  ### Can't do kmeans with a factor
  expect_equal(res[["splitCluster_twosplitsF_TRUE"]], NA)

  ## 10:         splitCluster       lvs2          TRUE
  ### Can't do kmeans with a factor
  expect_equal(res[["splitCluster_lvs2_TRUE"]], NA)

  ## 14:         splitCluster     constv          TRUE
  ### Can't use stop_splitby_constant=TRUE with a constant splitby
  expect_equal(res[["splitCluster_constv_TRUE"]], NA)

  ## 17:         splitCluster  twosplits         FALSE
  ### the splitting becomes random so results should differ between runs
  tabres1 <- with(res[["splitCluster_twosplits_FALSE"]], table(fin_grp, twosplits, exclude = c()))

  blah <- test_splitters_fn(
    sfn = "splitCluster", splitby = "twosplits",
    stopsplitting = FALSE
  )
  blah_det <- report_detections(blah, blockid = "bF")
  tabres2 <- with(blah_det, table(fin_grp, twosplits, exclude = c()))

  expect_false(isTRUE(all.equal(tabres1, tabres2)))
  ## 20:         splitCluster twosplitsF         FALSE
  ### Factors not allowed with splitCluster
  expect_equal(res[["splitCluster_twosplitsF_FALSE"]], NA)

  ## 24:         splitCluster       lvs2         FALSE
  expect_equal(res[["splitCluster_lvs2_FALSE"]], NA)

  ## 28:         splitCluster     constv         FALSE
  ### This next is random splits.
  tabres28 <- with(res[["splitCluster_constv_FALSE"]], table(fin_grp, constv, exclude = c()))

  ##  7: splitSpecifiedFactorMulti twosplitsF          TRUE
  ## splitSpecifiedFactorMulti requires a specific format for the labels of the factor
  expect_equal(res[["splitSpecifiedFactorMulti_twosplitsF_TRUE"]], NA)

  ## 11: splitSpecifiedFactor       lvs2          TRUE
  ### Only one non-zero entry in the following table per column because
  ### the blocks with specified factor level have to be kept together.
  resobj <- "splitSpecifiedFactor_lvs2_TRUE"
  tabres11 <- with(res[[resobj]], table(fin_grp, lvs2, exclude = c()))
  expect_true(all(apply(tabres11, 2, function(x) {
    sum(x > 0)
  })))
  expect_equal(uniqueN(res[[resobj]]$biggrp), uniqueN(res[[resobj]]$lvs2))

  ## 21: splitSpecifiedFactor twosplitsF         FALSE
  expect_equal(res[["splitSpecifiedFactor_twosplitsF_FALSE"]], NA)

  ## 25: splitSpecifiedFactor       lvs2         FALSE
  ### Starts within factor be then splits randomly
  tabres25 <- with(res[["splitSpecifiedFactor_lvs2_FALSE"]], table(fin_grp, lvs2, exclude = c()))

  theres2 <- findBlocks(
    idat = idat3, bdat = bdat4, blockid = "bF",
    pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla = Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE,
    splitfn = splitSpecifiedFactor, splitby = "lvs2", stop_splitby_constant = TRUE
  )
  theres2_det <- report_detections(theres2, blockid = "bF")
  expect_equal(uniqueN(theres2$biggrp), uniqueN(bdat4$lvs2))
  table(theres2$lvs2, theres2$biggrp)
})



##########################################################
## Testing node id functions here for now
nodeidfn <- function(d) {
  crc32hash <- digest::getVDigest(algo = "crc32")
  crc32hash(d)
}

i <- 1
node_id1 <- list()
node_id2 <- list()
for (i in 1:3) {
  if (i == 1) {
    node_id1[[i]] <- c(i * 2, i * 2 + 1)
    node_id2[[i]] <- nodeidfn(c(i, i + 1))
  } else {
    node_id1[[i]] <- lapply(node_id1[[i - 1]], function(j) {
      c(j * 2, j * 2 + 1)
    })
    node_id2[[i]] <- lapply(node_id2[[i - 1]], function(j) {
      nodeidfn(paste0(j, c(0, 1)))
    })
  }
}
