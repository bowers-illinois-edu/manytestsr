# Test and develop functions to control FDR via splitting
# context("FDR control for sequential structured testing")

library(data.table)
setDTthreads(1)
library(testthat)
library(here)
library(coin)
library(devtools)
library(onlineFDR)
load_all()
source(here::here("tests/testthat", "make_test_data.R"))

# Setup data

## This next is just for names
dists <- with(idat[bF == 1], dists_and_trans(Ynormb, Z))
## Blockwise tests
idat[, pvalue(independence_test(Ynormb ~ Z, data = .SD, teststat = "quadratic")), by = bF]

## Make some blocks have nothing going on
set.seed(12345)
idat[bF %in% c(7, 8), Ynormb := Ynull + rnorm(.N)]
idat[, pvalue(independence_test(Ynormb ~ Z, data = .SD, teststat = "quadratic")), by = bF]

idat[, c(names(dists), "rankY") := c(
  dists_and_trans(Ynormb, Z),
  list(Rfast::Rank(Ynormb)) # ,
), by = bF]

b1[, c(names(dists), "rankY") := c(
  dists_and_trans(Ynormb, Z),
  list(Rfast::Rank(Ynormb)) # ,
)]

Yvers <- c("Ynormb", names(dists), "rankY")
bigfmla_txt <- paste(paste(Yvers, collapse = "+"), "~Z|bF", sep = "")
bigfmla <- as.formula(bigfmla_txt)
blockfmla <- as.formula(paste(paste(Yvers, collapse = "+"), "~Z", sep = ""))

## Overall test
it1 <- independence_test(bigfmla, data = idat, teststat = "quadratic")
pvalue(it1)

## Blockwise tests
## Procedures below should try not to reject the null of no effects for three of the blocks.
pbs <- idat[, .(p = pvalue(independence_test(blockfmla, data = .SD, teststat = "quadratic"))), by = bF]
pbs$p

pbs_t <- idat[, .(p = pvalue(oneway_test(Ynormb ~ ZF, data = .SD, teststat = "quadratic"))), by = bF]
cbind(pbs, p_t = pbs_t$p)
cbind(pbs, p_t = pbs_t$p, sign(pbs$p - pbs_t$p)) ## interesting that t test is more powerful most of the time

## Using pre-specified splits. Set up the splits
## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
setkey(bdat4, bF)
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]
with(bdat4, table(lv1, lv3))
with(bdat4, table(lv1, lv2))
table(bdat4$lvs, exclude = c())

## This next allows a last split so use it for testing
bdat4[, gf1 := splitSpecifiedFactor(bF, x = lvs)]
with(bdat4, table(lvs, gf1))
bdat4[, gf2a := splitSpecifiedFactor(bF, x = lvs), by = gf1]
with(bdat4, table(lvs, gf2a, exclude = c()))
bdat4[, gf2 := interaction(gf1, gf2a)]
bdat4[, gf3a := splitSpecifiedFactor(bF, x = lvs), by = gf2]
with(bdat4, table(lvs, gf3a, exclude = c()))
bdat4[, gf3 := interaction(gf2, gf3a)]

## Does this look right? (yes for now)
bdat4[, .(lv1, lv2, lv3, lvs, gf1, gf2, gf3)]

## So now create a sequence of p-values
## Start with the global group (i.e. the group of all)
bdat4$gf0 <- factor("Z")

bdat4grps <- bdat4[, .(bF, gf0, gf1, gf2, gf3)]
idat_grps <- merge(idat, bdat4grps)
ftable(gf1 = idat_grps$gf1, gf2 = idat_grps$gf2, gf3 = idat_grps$gf3, exclude = c())

## This next just to simplify code
## Doesn't require monotonically increasing p-values here
simp_pfn <- function(fmla, dat) {
  pvalue(independence_test(fmla, data = dat, teststat = "quadratic"))
}

## Get the pvals for each split
for (i in 0:3) {
  gnm <- paste0("gf", i)
  idat_grps[, biggrp := get(gnm)]
  bdat4grps[, biggrp := get(gnm)]
  pb <- idat_grps[!is.na(biggrp), list(p = simp_pfn(fmla = bigfmla, dat = .SD)), by = biggrp]
  pnm <- paste0("p", i)
  setnames(pb, "p", pnm)
  setkeyv(pb, "biggrp")
  ## this next adds an i to the p-value so fix this
  ## bdat4grps[pb,,on="biggrp"]
  bdat4grps[pb, (pnm) := get(paste0("i.", pnm)), on = "biggrp"]
}

bdat4grps
pbs

bdat5 <- pbs[bdat4grps]

# Alpha Investing -- Foster 2008 style and/or Ramdas etc 2018 style
p_ex <- data.frame(
  id = 1:7,
  pval = c(.0001, .001, .004, .01, .1, .06, .04),
  group = c(
    "0", "0.1", "0.0", "0.1.0", "0.1.1",
    "0.0.0", "0.0.1"
  ),
  dateN = c(1, 2, 2, 3, 3, 3, 3)
)
p_ex$date <- as.character(p_ex$dateN)

## Alpha_investing doesn't need the whole set of hypotheses
res1 <- Alpha_investing(d = p_ex, date.format = "%d", random = FALSE)
res2 <- Alpha_investing(d = p_ex[1:3, ], date.format = "%d", random = FALSE)
res2a <- Alpha_investing(d = p_ex[1:3, ], date.format = "%d", random = TRUE)
## Date is just for batches or sequent, time in between dates is not relevant
res3 <- Alpha_investing(d = p_ex[1:3, ], date.format = "%Y", random = FALSE)
res4 <- Alpha_investing(d = p_ex[1:3, 1:2], random = FALSE) ## looks like the batching is not very relevant either although test 0.1 and 0.0 occur at same time
res5 <- Alpha_investing(d = p_ex[, 1:2], random = FALSE) ## batching becomes relevant in longer sequences
res6 <- Alpha_investing(d = p_ex[1, 1:2], random = FALSE) ## batching becomes relevant in longer sequences
res2a <- Alpha_investing(d = p_ex[c(1, 2, 4:5), ], date.format = "%d", random = FALSE)
res2b <- Alpha_investing(d = p_ex[c(1, 3, 6:7), ], date.format = "%d", random = FALSE)

## Effect of choosing starting alpha w_0
res1SAFFRON <- SAFFRON(d = p_ex, date.format = "%d", random = FALSE)
res2SAFFRON <- SAFFRON(d = p_ex[1:3, ], date.format = "%d", random = FALSE)
res6SAFFRON <- SAFFRON(d = p_ex[1, 1:2], random = FALSE) ## batching becomes relevant in longer sequences but the first alpha is set

res2SAFFRONa <- SAFFRON(d = p_ex[1:3, ], date.format = "%d", random = FALSE, w0 = .049)
res2a <- Alpha_investing(d = p_ex[1:3, ], date.format = "%d", random = FALSE, w0 = .049)



tmpdat <- data.frame(
  id = 8:15, pval = NA,
  group = paste(rep(p_ex$group[4:7], each = 2), rep(c(".0", ".1"), 4), sep = ""),
  dateN = rep(4:7, each = 2), date = NA
)
tmpdat$date <- as.character(tmpdat$dateN)
tmpdat$pval <- c(.02, .03, .2, .01, .07, .01, .05, .01)

p_ex2 <- rbind(p_ex, tmpdat)



# p0 root
a0 <- Alpha_investing(d = p_ex2[1, ], date.format = "%d", random = FALSE, w0 = .025)
# next level, groups under group 0:  0.1 and 0.0
a1 <- Alpha_investing(d = p_ex2[c(1, 2:3), ], date.format = "%d", random = FALSE)
## then groups under 0.1 : 0.1.0 and 0.1.1
a21 <- Alpha_investing(d = p_ex2[c(1, 2, 4:5), ], date.format = "%d", random = FALSE)
## groups under 0.1.0: 0.1.0.0 and 0.1.0.1
a210 <- Alpha_investing(d = p_ex2[c(1, 2, 4, 8:9), ], date.format = "%d", random = FALSE)
## groups under 0.1.1: 0.1.1.0 and 0.1.1.1
a211 <- Alpha_investing(d = p_ex2[c(1, 2, 5, 10:11), ], date.format = "%d", random = FALSE)
## groups under 0.0: 0.0.0 and 0.0.1
a20 <- Alpha_investing(d = p_ex2[c(1, 3, 6:7), ], date.format = "%d", random = FALSE)
## groups under 0.0.0: 0.0.0.0 and 0.0.0.1
a200 <- Alpha_investing(d = p_ex2[c(1, 3, 6, 12:13), ], date.format = "%d", random = FALSE)
## groups under 0.0.1: 0.0.1.0 and 0.0.1.1
a201 <- Alpha_investing(d = p_ex2[c(1, 3, 7, 14:15), ], date.format = "%d", random = FALSE)

## library(debug)
load_all()
debug(findBlocks)
blah <- findBlocks(
  idat = idat, bdat = bdat4, blockid = "bF", splitfn = splitSpecifiedFactor,
  pfn = pIndepDist, alphafn = alpha_investing, simthresh = 20,
  sims = 1000, maxtest = 100, thealpha = 0.05,
  fmla = Ynormb ~ ZF | bF,
  parallel = "multicore", copydts = TRUE, splitby = "lvs"
)


break

library(data.table)
library(stringi)
load("temp_pb_tracker.rda")
pbtracker$alpha4 <- NULL
pbtracker$biggrpC_prev <- NULL
pbtracker[, depth := as.numeric(stri_sub(batch, 2, 2))]


## Can we vectorize?
lv <- 4
mid_roots <- pbtracker[depth == (lv - 1) & (testable), nodenum]
# paths <- matrix(NA, nrow = (lv + 1), ncol = length(mid_roots))
## Each mid root has two children
# paths[(lv - 1):(lv + 1), ] <- vapply(mid_roots, FUN.VALUE=c(0,0,0),FUN=function(x) {
# c(x, 2 * x, 2 * x + 1)
# })
paths <- lapply(mid_roots, lv = lv, FUN = function(x, lv) {
  ## Each mid root has two children
  kids <- c(x, 2 * x, 2 * x + 1) ## mid_root is part of kids vector
  ## But each mid root may have many parents
  parents <- rep(NA, length = lv - 2)
  parents[lv - 2] <- floor(x / 2) # l = level - 1
  l <- lv - 3
  while (l > 0) {
    parents[l] <- floor(min(parents, na.rm = TRUE) / 2)
    l <- l - 1
  }
  return(c(parents, kids))
})

## But each mid root may have many parents
## for (l in seq((lv - 1), 2)) {
##  paths[l - 1, ] <- floor(paths[l, ] / 2)
## }
##

pbtracker[10, p := .2] ## testing the function
setkeyv(pbtracker, "nodenum")
for (j in 1:(lv - 1)) {
  pbtracker[paths[j], alpha4 := alpha_investing(pval = p, batch = batch)]
}

pbtracker[paths[1], alpha_investing(pval = p, batch = batch)]
pbtracker[paths[2], alpha_investing(pval = p, batch = batch)]


break
#####
pbtracker[nodenum %in% paths2[, 1]]
pbtracker[nodenum %in% paths2[, 2]]
pbtracker[nodenum %in% paths2[, 3]]


nodes_lv <- pbtracker[nodenum <= 2^lv & nodenum >= 2^(lv - 1), nodenum]
paths <- matrix(NA, nrow = lv, ncol = length(nodes_lv))
paths[lv, ] <- nodes_lv
for (l in seq(lv, 2)) {
  paths[l - 1, ] <- floor(paths[l, ] / 2)
}
roots <- unique(paths[1:(lv - 1), ], MARGIN = 2)







pb_lv_4 <- pbtracker[nodenum <= 2^lv & nodenum >= 2^(lv - 1), ]
pb_lv_4[, par_level_4 := floor(nodenum / 2)]
pb_lv_4[, par_level_3 := floor(par_level_4 / 2)]
pb_lv_4[, par_level_2 := floor(par_level_3 / 2)]
## Still need lists, I think.

pb_lv_4[par_level_4 == 4, .(unique(c(nodenum, par_level_4, par_level_3, par_level_2)))]



find_path <- function(nodenum) {
  path <- vector(mode = "integer")
  parent_node <- nodenum
  path[1] <- parent_node
  while (parent_node > 1) {
    parent_node <- floor(parent_node / 2)
    path <- c(path, parent_node)
  }
  return(sort(path))
}
###
tree_level <- 1
pbtracker[nodenum == tree_level, ]
tree_level <- 2
pbtracker[nodenum %in% find_path(nodenum[3]), ]

tree_level <- 4
the_leaves <- pbtracker[nchar(biggrpC) == tree_level, ]

tree_level_prev <- tree_level - 1
roots_current <- pbtracker[nchar(biggrpC) == tree_level_prev, ]

testing_paths_4 <- lapply(roots_current$nodenum, function(nodenum) {
  partial_path <- c(1, floor(nodenum / 2), nodenum, c(2 * nodenum, 2 * nodenum + 1))
  return(partial_path)
})

testing_paths_5 <- lapply(the_leaves$nodenum, function(nodenum) {
  partial_path <- c(1, floor(nodenum / 2), nodenum, c(2 * nodenum, 2 * nodenum + 1))
  return(partial_path)
})


pbtracker[str_sub(biggrpC, 1, 2) == biggrpC[3], ]
pbtracker[str_sub(biggrpC, 1, 3) == biggrpC[4], ]

g <- pbtracker$biggrpC[2]
## How far can we go with g? And then the next one?

find_path_from_leaf <- function(leaf_string) {
  leaf_depth <- nchar(leaf_string)
  path_string <- vector(mode = "character", length = leaf_depth)
  for (i in 1:leaf_depth) {
    path_string[i] <- stri_sub(leaf_string, 1, i)
  }
  return(path_string)
}

path_to_8 <- find_path_from_leaf(pbtracker$biggrpC[8])
path_to_9 <- find_path_from_leaf(pbtracker$biggrpC[9])


children_of_4 <- c(2 * 4, 2 * 4 + 1)
pbtracker[pbtracker$nodenum %in% children_of_4, ]
parent_of_4_and_5 <- unique(floor(c(4, 5) / 2))

# FYI children are "index*2" and "index*2+1
## For any node x -> i It's children will be 2*i and 2*i + 1 And it's parent will be floor(i/2)
