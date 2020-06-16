# Test and develop functions to split the data
context("Performance of Splitting Functions")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- TRUE
if(interactive){
 library(here)
 library(data.table)
 source(here::here("tests/testthat", "make_test_data.R"))
 devtools::load_all() ## use  this during debugging
}

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
## for the first split (say, with only two groups or even with only one group), splitLOO will choose one block at random. So, splitLOO with a splitting vector/criteria/variable that does not vary is the same as random splits. So, probably ok for error rate but not great for power.

## A splitting variable with little variation.
### START HERE
set.seed(12345)
bdat4[,twosplits:=rbinom(.N,1,.5)]
bdat4[,twosplitsF:=factor(twosplits)] ## should only have 2 splits
bdat4[,lvs2:=interaction(lv1,lv2)] ## should only have 4 spits
bdat4[,constv:=rep(10,.N)]

 theres1 <- findBlocks(idat = idat3, bdat = bdat4, blockid = "bF",     pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla =  Ytauv2 ~ ZF | bF,
    parallel = 'no', copydts = TRUE,
    splitfn = splitCluster, splitby = 'twosplits', stop_splitby_constant = TRUE
  )
theres1_det <- report_detections(theres1,blockid="bF")
expect_equal(uniqueN(theres1$biggrp),uniqueN(bdat4$twosplits))
## This next because I know that we should reject in both groups.
expect_equal(uniqueN(theres1_det$fin_grp),uniqueN(bdat4$twosplits))
table(theres1$twosplitsF,theres1$biggrp)

## Not testing "splitSpecified" because splitSpecifiedFactor does a better job.
split_test_parms <- data.table(expand.grid(sfn = c("splitLOO", "splitEqualApprox",
                                         "splitCluster", "splitSpecifiedFactor"),
                                splitby = c("twosplits","twosplitsF","lvs2","constv"),
                                stopsplitting =  c(TRUE,FALSE), stringsAsFactors=FALSE))

split_test_parms <- split_test_parms[sfn!="splitSpecifiedFactor" | 
                                     (sfn=="splitSpecifiedFactor" & !(splitby %in% c("twosplits","constv"))),]

test_splitters_fn <- function(sfn,splitby,stopsplitting){
    theres <- findBlocks(idat = idat3, bdat = bdat4, blockid = "bF",
                         pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
                         fmla =  Ytauv2 ~ ZF | bF,
                         parallel = 'no', copydts = TRUE,
                         splitfn = get(sfn), splitby = splitby, stop_splitby_constant = stopsplitting
    )
    return(theres)
}

res <- mapply(FUN=function(sfn,sby,stopsplitting){
                  message(paste(sfn, sby, stopsplitting, collapse = ","))
                  try(test_splitters_fn(sfn=sfn,
                                    splitby=sby,
                                    stopsplitting=stopsplitting))
    },
    sfn=split_test_parms$sfn,
    sby=split_test_parms$splitby,
    stopsplit=split_test_parms$stopsplitting,
    SIMPLIFY=FALSE)


## splitLOO will tend to first choose one of twosplits==1 at random and then again until there are no more twosplits==1, then it will stop. If it would reject with the entire group of twosplits==0, it will still stop but then report that as a rejection.
blah1 <- test_splitters_fn(sfn="splitLOO",splitby="twosplits",stopsplitting=TRUE)
blah1_det <- report_detections(blah1,blockid="bF")
blah1_tree <- make_graph(make_tree(blah1))

blah2 <- test_splitters_fn(sfn="splitLOO",splitby="twosplits",stopsplitting=FALSE)
blah2_det <- report_detections(blah2,blockid="bF")
blah2_tree <- make_graph(make_tree(blah2))

## This next should stop after one split --- returning two groups that map directly onto twosplits
blah3 <- test_splitters_fn(sfn="splitCluster",splitby="twosplits",stopsplitting=TRUE)
blah3_det <- report_detections(blah3,blockid="bF")
blah3_tree <- make_graph(make_tree(blah3))

## This should continue splitting although I think that all of the splitting should still be within twosplits
blah4 <- test_splitters_fn(sfn="splitCluster",splitby="twosplits",stopsplitting=FALSE)
blah4_det <- report_detections(blah4,blockid="bF")
blah4_tree <- make_graph(make_tree(blah4))

## This next should stop after 4 splits.
blah5 <- test_splitters_fn(sfn="splitSpecifiedFactor",splitby="lvs2",stopsplitting=TRUE)
blah5_det <- report_detections(blah5,blockid="bF")
blah5_tree <- make_graph(make_tree(blah5))
with(blah5_det,table(hit,lvs2,exclude=c()))
with(blah5_det,table(fin_grp,lvs2,exclude=c()))

## This should continue splitting although I think that all of the splitting should still be within lvs2
blah6 <- test_splitters_fn(sfn="splitSpecifiedFactor",splitby="lvs2",stopsplitting=FALSE)
blah6_det <- report_detections(blah6,blockid="bF")
blah6_tree <- make_graph(make_tree(blah6))
with(blah6_det,table(hit,lvs2,exclude=c()))
with(blah6_det,table(fin_grp,lvs2,exclude=c()))


##test_that("Splitters stop when given a splitting criteria vector with few values and stopping is requested.", {
##
 ##   })





theres2 <- findBlocks(idat = idat3, bdat = bdat4, blockid = "bF",
                      pfn = pIndepDist, alphafn = NULL, thealpha = 0.05,
    fmla =  Ytauv2 ~ ZF | bF,
    parallel = 'no', copydts = TRUE,
    splitfn = splitSpecifiedFactor, splitby = 'lvs2', stop_splitby_constant = TRUE
  )
theres2_det <- report_detections(theres2,blockid="bF")
expect_equal(uniqueN(theres2$biggrp),uniqueN(bdat4$lvs2))
table(theres2$lvs2,theres2$biggrp)

### START HERE OR JUST ABOVE TO VERIFY THAT THIS WORKS FOR ALL SPLITTERS. ALSO CHECK FOR BEHAVIOR WHEN stop_splitby_constant=FALSE

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
