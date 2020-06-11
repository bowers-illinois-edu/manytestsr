# Test and develop functions to split the data
context("Performance of Splitting Functions")

# library(here)
# library(data.table)
# source(here::here("tests/testthat", "make_test_data.R"))
# devtools::load_all() ## use  this during debugging

## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]

## The test of leave one out splitting is easiest.
testSplitLOO <- splitLOO(bdat4$bF, bdat4$hwt)
stopifnot(sum(testSplitLOO == "1") == length(testSplitLOO) - 1)

## Not sure what kind of test to write here. Leaving this as is for now. Hoping it doesn't break package building
bdat4[, g1 := splitCluster(bid = as.character(bF), x = hwt)]
bdat4[, g2 := splitEqualApprox(bid = as.character(bF), x = hwt)]
bdat4[, g3 := splitEqual(bid = as.character(bF), x = hwt)]
bdat4[, g4 := splitEqual(bid = as.character(bF), x = v1)]
bdat4[, g5 := splitCluster(bid = as.character(bF), x = v4)]

# Maybe we want something from splitCLuster like the difference between the  means is bigger than any random split?
bdat4[, .(mn = mean(hwt), sd = sd(hwt)), by = g1]
bdat4[, .(mn = mean(hwt), sd = sd(hwt)), by = g2]
bdat4[, .(mn = mean(hwt), sd = sd(hwt)), by = g3]
bdat4[, .(mn = mean(v1), sd = sd(v1)), by = g4]

## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]

with(bdat4, table(lv1, lv3))
with(bdat4, table(lv1, lv2))
table(bdat4$lv1_2)

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

set.seed(12355)
idat$y1test_null <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 0, prop_blocks_0 = .5)
idat$y1test_zeros <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 2, prop_blocks_0 = .5)
idat[, Y_zeros := y1test_zeros * Z + y0 * (1 - Z)]
idat[, Y_null := y1test_null * Z + y0 * (1 - Z)]


#### When splitby has no variation within biggrp, further splitting doesn't do anything and the algorithm may just continue to produce the same p-values until maxtest is reached depending on the splitting function. So, for example, in splitSpecifiedFactor and splitCluster we want the algorithm to stop or switch to another splitter (TODO maybe) once clusters / branches have no variation on splitby  --- any further splitting would be essentially random.

## In splitSpecifiedFactor we want it to stop when splitby is constant within group
## In splitCluster we want it to stop when splitby is constant within group or move to random splits or splitLOO or splitEqualApprox
## In splitEqualApprox we want splitby to be something like N and so we want it to continue even if the groups have the same N.
## In splitLOO could stop when splitby is constant within the bigger group or start choosing blocks to focus on at random.

## A splitting variable with little variation.
### START HERE
set.seed(12345)
bdat4[,twosplits:=rbinom(.N,1,.5)]
bdat4[,twosplitsF:=factor(twosplits)]

splittingfns <- c("splitLOO", "splitEqualApprox", "splitCluster", "splitSpecifiedFactor")

 theres1 <- findBlocks(idat = idat3, bdat = bdat4, blockid = "bF", splitfn = splitCluster,
    pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla =  Ytauv2 ~ ZF | bF,
    parallel = "no", copydts = TRUE, splitby = 'twosplitsF'
  )


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
