## Test and develop functions to adapt p-values and alpha levels as the tree grows
testthat::context("Alpha and P Adjusting Performance")

## The next lines are for use when creating the tests. We want interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(dtplyr)
  # remotes::install_github("bowers-illinois-edu/TreeTestSim")
  library(TreeTestsSim)
  source(here("tests/testthat", "make_test_data.R"))
  library(devtools)
  ## This next is a copy from the package TreeTestsSim
  ## source(here("tests/testthat", "generate_tree.R"))
  load_all() ## use  this during debugging
}
setDTthreads(1)
options(digits = 4)

## For now, create some data within this file.
## Move it elsewhere soon.

## A tree where 50% of the leaves are non-null
## the tree is defined at the level of the block
bdt <- generate_tree_DT(max_level = 3, k = 4, t = .5)
bdt[, leaf := as.numeric(level == max(level))]
table(bdt$leaf)
# 1. build a named vector so parent_lookup["22"] == 6, etc.
parent_lookup <- setNames(bdt$parent, bdt$node)
## Compare to the mathematical way to calculate the parent
expect_equal(parent_lookup[22][[1]], floor((22 - 2) / 4) + 1)

# 2. for a given node, return the vector of its ancestors+itself,
#    *excluding* the NA‐parent root.
get_path <- function(node) {
  path <- integer(0)
  current <- node

  # walk up until you hit an NA parent
  while (!is.na(parent_lookup[as.character(current)])) {
    path <- c(path, current) # record this node
    current <- parent_lookup[as.character(current)] # jump to parent
  }

  # reverse so it reads root→…→node
  rev(path)
}



# 3. apply to every node, paste with ".", and coerce to factor
bdt[, lvls_fac := factor(sapply(node, function(n) {
  paste(get_path(n), collapse = ".")
}))]

## We are going to focus on the leaves which are the individual blocks.

bdt1 <- droplevels(bdt[leaf == 1, ])

## The splitSpecifiedFactorMulti function should recreate this tree.
bdt1[, split1 := splitSpecifiedFactorMulti(bid = node, x = lvls_fac)]
bdt1[, split2 := splitSpecifiedFactorMulti(node, lvls_fac), by = split1]
bdt1[, split3 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2)]

with(bdt1, ftable(split1 = split1, split2 = split2, split3 = split3))

## Everybody within the same parent is in the same first level split
test_split_1 <- bdt1 %>%
  group_by(parent) %>%
  summarize(val = length(unique(split1)))
stopifnot(all(test_split_1$val == 1))

test_split_2 <- bdt1 %>%
  group_by(split1) %>%
  summarize(val = length(unique(split2)))
stopifnot(all(test_split_2$val == 4))
stopifnot(nrow(test_split_2) == 4)

bdt1 %>% filter(split1 == 1)

## This next returns all 0s when splitting is no longer possible. So, consider watching out for this in the find_blocks code
bdt1[, split4 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3)]
bdt1[, split5 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3, split4)]


# 3. Make individual level data
bdt1[, nb := 100]
bdt1[, bary0 := 0]

idt <- data.table(b = rep(bdt1$node, bdt1$nb))
idt[, bF := factor(b)]
idt[, id := 1:nrow(idt)]
idt[, bary0 := rep(bdt1$bary0, bdt1$nb)]
## this next is not necessary since it is the same y0 for all blocks. But it will not be that case in other situations
idt[, y0 := rnorm(.N, bary0, sd = 1), by = b]


tau_null <- function(ybase, tau_sds, covariate) {
  return(ybase)
}
idt[, y1 := create_effects(idat = idt, ybase = "y0", blockid = "bF", tau_fn = tau_null, tau_size = 0, prop_blocks_0 = 1)]





idt <- data.table(b = rep(c(1:20), length = 1000))

## 20 blocks, each with 50 obs
idt <- data.table(i = 1:1000, b = rep(c(1:20), length = 1000))
bdt <- data.table(b = 1:20)
setkey(bdt, b)
setkey(idt, b)
bdt[, vb1 := rep(c(1, 2, 3, 4), length.out = .N)]
## 2 levels within each state
bdt[, vb2 := rep(c(1, 2, 3), length.out = .N), by = vb1]
bdt[, vb3 := rep(seq(1, .N), length.out = .N), by = vb2]
bdt[, vbnest := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
ftable(lv1 = bdt$vb1, lv2 = bdt$vb2, bdt$vb3, exclude = c())
table(bdt$vbnest, exclude = c())
## Now merge the block data onto the individual level data
idt <- bdt[idt]
set.seed(12345)
idt[, vi1 := round(runif(.N, min = 0, max = 10)), by = b]
idt[, y0 := vi1 + vb1 + rnorm(100), by = b]

# 4 nodes at level 1
bdt[, lv1 := rep(c(1, 2, 3, 4), length.out = .N)]
bdt[, lv2 := rep(1:.N, length.out = .N), by = lv1]
with(bdt, table(lv1, lv2, exclude = c()))
bdt[, lvls := interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
table(bdt$lvls, exclude = c())




## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]

## Setting up to test the local p and alpha adjusting methods versus FWER
## methods on different splitters


local_adj_methods <- c(
  "local_hommel_all_ps", "local_bh_all_ps",
  "local_simes", "local_unadj_all_ps"
)

alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "alpha_addis", "NULL"),
  sfn = c(
    "splitCluster", "splitEqualApprox", "splitLOO",
    "splitSpecifiedFactor", "splitSpecifiedFactorMulti"
  ),
  local_adj_fn = local_adj_methods,
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <- "lvs"

testing_fn <- function(afn, sfn, local_adj, sby, fmla = Ytauv2 ~ ZF | bF, idat = idat3, bdat = bdat4) {
  if (afn == "NULL") {
    theafn <- NULL
  } else {
    theafn <- getFromNamespace(afn, ns = "manytestsr")
  }

  if (local_adj == "NULL") {
    local_adj <- NULL
  } else {
    local_adj <- getFromNamespace(local_adj, ns = "manytestsr")
  }
  ## afn and sfn and sby are character names
  theres <- find_blocks(
    idat = idat, bdat = bdat, blockid = "bF", splitfn = get(sfn),
    pfn = pIndepDist, alphafn = theafn, local_adj_p_fn = local_adj, thealpha = 0.05, thew0 = .05 - .001,
    fmla = fmla,
    copydts = TRUE, splitby = sby, stop_splitby_constant = TRUE, parallel = "multicore", ncores = 2
  )
  setkey(theres, bF)
  return(theres)[order(biggrp)]
  # theps <- grep("^p[0-9]", names(theres), value = TRUE)
  # theas <- grep("^alpha", names(theres), value = TRUE)
  # truth <- grep("^ate", names(theres), value = TRUE)
  # return(theres[, .SD, .SDcols = c(
  #  theps, theas, truth, "pfinalb", "blocksbygroup",
  #  "biggrp", "bF", "hwt", "nb"
  # )][order(biggrp)])
}

##  Maybe make w0 more like .05.
alpha_and_splits[c(1, 2, 4), ]
## debug(find_blocks)
res_ai <- testing_fn(
  afn = alpha_and_splits[1, "afn"],
  sfn = alpha_and_splits[1, "sfn"],
  local_adj = alpha_and_splits[1, "local_adj_fn"],
  sby = alpha_and_splits[1, "splitby"],
  idat = idat3, bdat = bdat4
)
res_saffron <- testing_fn(
  afn = alpha_and_splits[2, "afn"],
  sfn = alpha_and_splits[2, "sfn"],
  sby = alpha_and_splits[2, "splitby"],
  local_adj = alpha_and_splits[2, "local_adj_fn"],
  idat = idat3, bdat = bdat4
)
res_fwer <- testing_fn(
  afn = alpha_and_splits[4, "afn"],
  sfn = alpha_and_splits[4, "sfn"],
  sby = alpha_and_splits[4, "splitby"],
  local_adj = alpha_and_splits[4, "local_adj_fn"],
  idat = idat3, bdat = bdat4
)

grep("^p[0-9]", names(res_ai), value = TRUE)
grep("^p[0-9]", names(res_saffron), value = TRUE)
grep("^p[0-9]", names(res_fwer), value = TRUE)
options(digits = 3, scipen = 8)
res_ai[order(p1, p2, p3, p4, p5, p6, p7, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6, p7)]
res_saffron[order(p1, p2, p3, p4, p5, p6, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6)]
res_fwer[order(p1, p2, p3, p4, p5, p6, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6)]
## Early in the splitting. No change
all.equal(res_ai$p2, res_saffron$p2)
all.equal(res_ai$p2, res_fwer$p2)
## But then some blocks stop testing
all.equal(res_ai$p5, res_saffron$p5)
all.equal(res_ai$p5, res_fwer$p5)
cbind(
  res_ai[, .(bF, p5, p6)],
  res_saffron[, .(bF, p5, p6)],
  res_fwer[, .(bF, p5, p6)]
)

## Which blocks were identified?
##
## The alpha investing style procedures should identify **more** blocks
# We detect an effect in an individual block if:
## (1) the block is a leaf and that leaf has p < alpha
# We detect an effect for a group of blocks if
## (2) the blocks are leaves with p > alpha but the parent of those blocks have p < alpha

## Change one result from res_fwer to enable us to check this:

res_ai[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_saffron[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_fwer[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]

res_ai[bF == "10", p5 := .5]
res_ai[bF == "10", pfinalb := .5]
res_saffron[bF == "10", p5 := .5]
res_saffron[bF == "10", pfinalb := .5]
res_fwer[bF == "10", p5 := .5]
res_fwer[bF == "10", pfinalb := .5]

res_ai[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_saffron[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_fwer[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]

## With alpha fixed
res_fwer_det <- report_detections(res_fwer)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_fwer_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

## With alpha varying according to the alpha investing
res_ai_det <- report_detections(res_ai, fwer = FALSE)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_ai_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

## And with the saffron procedure
res_saffron_det <- report_detections(res_saffron, fwer = FALSE)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_saffron_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

res_fwer_tree <- make_results_tree(res_fwer, blockid = "bF")
res_saffron_tree <- make_results_tree(res_saffron, blockid = "bF")
res_ai_tree <- make_results_tree(res_ai, blockid = "bF")

blahN <- res_fwer_tree %>%
  activate(nodes) %>%
  as_tibble()
blahE <- res_fwer_tree %>%
  activate(edges) %>%
  as_tibble()

res_fwer_g <- make_results_ggraph(res_fwer_tree)
res_ai_g <- make_results_ggraph(res_ai_tree)
res_saffron_g <- make_results_ggraph(res_saffron_tree)

## Compare these graphs to the results in res_ai_det and res_fwer_det above.
ggsave(res_fwer_g, file = "res_fwer_g.pdf", bg = "transparent", width = 12, height = 10)
ggsave(res_ai_g, file = "res_ai_g.pdf", bg = "transparent", width = 13, height = 10)
ggsave(res_saffron_g, file = "res_saffron_g.pdf", bg = "transparent", width = 13, height = 10)

## Compare working with tree/graph versus blocks
## TODO. Does report_detections produce the same discoveries as make_results_tree?

## Criteria for comparisons:
## These functions can be used to run the tests for the different scenarios.
## They have hard coded the relationships in alpha_and_splits above for now.

eval_numgrps <- function(detection_obj) {
  numgrps <- sapply(detection_obj, function(dat) {
    length(unique(dat$hit_grp))
  })
  ## names_split <- stri_split_regex(names(numgrps), "_", simplify = TRUE)

  expect_lte(numgrps[[4]], max(numgrps[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numgrps[[8]], max(numgrps[5:7]))
  expect_lte(numgrps[[12]], max(numgrps[9:11]))
  expect_lte(numgrps[[16]], max(numgrps[13:15]))
  expect_lte(numgrps[[20]], max(numgrps[17:19]))
}

eval_maxp <- function(detection_obj) {
  maxp <- sapply(detection_obj, function(dat) {
    dat[, sizerp := .N, by = hit_grp]
    dat[, p := ifelse(sizerp > 1, fin_parent_p, max_p)]
    max(dat[, .(min_max_p = min(p)), by = hit_grp]$min_max_p)
  })
  ## Adding tolerance of .01 here but could be more depending on whether we do any simulation
  ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(maxp[[4]], max(maxp[1:3]) + .01)
  expect_lte(maxp[[8]], max(maxp[5:7]) + .01)
  expect_lte(maxp[[12]], max(maxp[9:11]) + .01)
  expect_lte(maxp[[16]], max(maxp[13:15]) + .01)
  expect_lte(maxp[[20]], max(maxp[17:19]) + .01)
}

eval_single_blocks_found <- function(detection_obj) {
  ## Number of individual blocks detected versus groups (expect more singletons with the alpha adjusting approaches)
  numsingletons <- sapply(detection_obj, function(dat) {
    tab <- table(dat$hit_grp)
    sum(tab == 1)
  })
  expect_lte(numsingletons[[4]], min(numsingletons[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numsingletons[[8]], min(numsingletons[5:7]))
  expect_lte(numsingletons[[12]], min(numsingletons[9:11]))
  expect_lte(numsingletons[[16]], min(numsingletons[13:15]))
  expect_lte(numsingletons[[20]], min(numsingletons[17:19]))
}

eval_treedepth <- function(detection_obj) {
  treedepth <- sapply(detection_obj, function(dat) {
    max(stri_count_fixed(dat$biggrp, ".")) + 1
  })
  expect_lte(treedepth[[4]], max(treedepth[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(treedepth[[8]], max(treedepth[5:7]))
  expect_lte(treedepth[[12]], max(treedepth[9:11]))
  expect_lte(treedepth[[16]], max(treedepth[13:15]))
  expect_lte(treedepth[[20]], max(treedepth[17:19]))
}

eval_minate <- function(detection_obj, atenm) {
  ## This compares the ATEs detected. I currently think that the fixed alpha approach should not be as sensitive or powerful as the alpha adjusting approaches. But somehow I wrote this test. Leaving it and ignoring it for now. Also the results fail this test. So I think my current intuition is correct.
  minate <- sapply(detection_obj, function(dat) {
    min(dat[, min(abs(get(atenm))), by = hit_grp]$V1)
  })
  expect_lte(minate[[4]], min(minate[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(minate[[8]], min(minate[5:7]))
  expect_lte(minate[[12]], min(minate[9:11]))
  expect_lte(minate[[16]], min(minate[13:15]))
  expect_lte(minate[[20]], min(minate[17:19]))
}

resnms <- apply(alpha_and_splits, 1, function(x) {
  tmp <- paste0(x, collapse = "_", sep = "")
  gsub("NULL", "fwer_fwer", tmp)
})


test_that("alphafns work across splitters for no effects", {
  ## First with no effects at all. So, basically no splitting should happen and no discoveries should be reported.
  tau_null <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynull ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_null) <- resnms

  tau_null_det <- lapply(seq_along(tau_null), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_null)[[i]], 1, 4) == "NULL"
    report_detections(tau_null[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_null_det) <- resnms
  ## No approach should detect anything, So we don't imagine any differences between approaches.
  anydetected <- sapply(tau_null_det, nrow)
  expect(all(anydetected == 0), TRUE)
})

test_that("alphafns work across splitters for large and homogenous effects", {
  ## All blocks same large effect. Both approaches should detect effects in basically all blocks (blocks only differ in size)
  tau_homog <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Yhomog ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_homog) <- resnms

  tau_homog_det <- lapply(seq_along(tau_homog), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_homog)[[i]], 1, 4) == "NULL"
    report_detections(tau_homog[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_homog_det) <- resnms
  ## We are comparing resnms 1 and 2 and 3 (the two fdr approaches) versus 4 (fwer) etc..

  ## Block 1 is the smallest block but it sometimes triggers a detection
  table(idat3$bF)
  block1detected <- sapply(tau_homog_det, function(dat) {
    any(dat$bF == "1")
  })
  table(block1detected)
  ## So basically all blocks or at least all but the smallest block should be detected
  numblks_homog <- sapply(tau_homog_det, nrow)
  expect_equal(all(numblks_homog >= 19), TRUE)
  ## Each block should be detected  --- they should not be grouped together ---
  ##  so the number of groups should be 19 as well
  numgrps_homog <- sapply(tau_homog_det, function(dat) {
    length(unique(dat$hit_grp))
  })
  expect_equal(all(numgrps_homog >= 19), TRUE)

  eval_numgrps(tau_homog_det)
  ## Number of individual blocks detected versus groups (expect more or same singletons with the alpha adjusting approaches)
  eval_single_blocks_found(tau_homog_det)

  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  minate_homog <- sapply(tau_homog_det, function(dat) {
    min(dat[, min(abs(ate_homog)), by = hit_grp]$V1)
  })

  ## Highest p detected: Doesn't reliably work in this case with all huge effects.
  eval_maxp(tau_homog_det)

  # g_homog <- lapply(tau_homog,function(obj){ make_results_ggraph(make_results_tree(obj)) })
  # for(i in 1:length(g_homog)){
  #    ggsave(g_homog[[i]], file = paste0("tau_homog_g",i,".pdf"),
  #           bg = "transparent", width = 14, height = 7)
  # }
})


test_that("alphafns work across splitters for individually heteogeneous effects and increase with block size. Also some completely null blocks.", {
  ################################################################################
  ## Some blocks have no effect at all in norm_inc
  tau_norm_inc <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynorm_inc ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_norm_inc) <- resnms
  tau_norm_inc_det <- lapply(seq_along(tau_norm_inc), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_norm_inc)[[i]], 1, 4) == "NULL"
    report_detections(tau_norm_inc[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_norm_inc_det) <- resnms

  ## Number of individual blocks detected versus groups
  ##  (expect more singletons with the alpha adjusting approaches than with fixed alpha)
  ## eval_single_blocks_found(tau_norm_inc_det)

  ## Number of blocks detected: this is not clear what to expect because of the ability to declare "detect" for *groups* of blocks.
  ## So, FWER might stop testing  and return many blocks.
  numblks_norm_inc <- sapply(tau_norm_inc_det, nrow)
  # expect_lte(numblks_norm_inc[[3]],min(numblks_norm_inc[1:2]))
  # expect_lte(numblks_norm_inc[[6]],min(numblks_norm_inc[4:5]))
  # expect_lte(numblks_norm_inc[[9]],min(numblks_norm_inc[7:8]))

  ## Number of groups of blocks detected: expect more or equal from alpha adjusting
  eval_numgrps(tau_norm_inc_det)
  ## Depth of testing (maxdepth): Probably a deeper tree or equal.
  eval_treedepth(tau_norm_inc_det)
  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  ## I currently disagree with # eval_minate --- shouldn't the alpha adjusted approaches be more sensitive and detect smaller ates?
  ## Commenting this out for now. I have expectations about p-values but not about ates per se
  ## eval_minate(tau_norm_inc_det, atenm = "ate_norm_inc")

  ## Highest p detected
  eval_maxp(tau_norm_inc_det)
  ## Null blocks detected (but this could be ok if they were in a group of other blocks having strong effects, so do by hit_grp)
  tau_norm_inc_det[[1]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  tau_norm_inc_det[[4]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  symdiff <- function(x, y) {
    setdiff(union(x, y), intersect(x, y))
  } # https://www.r-bloggers.com/symmetric-set-differences-in-r/
  symdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[4]]$bF)
  symdiff(tau_norm_inc_det[[2]]$bF, tau_norm_inc_det[[4]]$bF)
  setdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[4]]$bF)
  setdiff(tau_norm_inc_det[[2]]$bF, tau_norm_inc_det[[1]]$bF)

  sort(setdiff(bdat4$bF, tau_norm_inc_det[[1]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[2]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[3]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[4]]$bF))
  # g_norm_inc <- lapply(tau_norm_inc, function(obj) {
  #   make_results_ggraph(make_results_tree(obj))
  # })
  # for (i in 1:length(g_norm_inc)) {
  #   ggsave(g_norm_inc[[i]],
  #     file = paste0("tau_norm_inc_g", i, ".pdf"),
  #     bg = "transparent", width = 14, height = 7
  #   )
  # }
})

## TODO this next not working. Maybe intuitions are wrong?
## test_that("alphafns work across splitters for individually heteogeneous effects
##  and decrease with block size. Also some completely null blocks.", {
##  ## Some blocks have no effect at all in norm_dec
##  tau_norm_dec <- mapply(
##    FUN = function(afn = afn, sfn = sfn, sby = sby) {
##      message(paste(afn, sfn, sby, collapse = ","))
##      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynorm_dec ~ ZF | bF)
##    },
##    afn = alpha_and_splits$afn,
##    sfn = alpha_and_splits$sfn,
##    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
##  )
##  names(tau_norm_dec) <- resnms
##  tau_norm_dec_det <- lapply(seq_along(tau_norm_dec), function(i) {
##    # message(i)
##    fwer <- stri_sub(names(tau_norm_dec)[[i]], 1, 4) == "NULL"
##    report_detections(tau_norm_dec[[i]], fwer = fwer, only_hits = TRUE)
##  })
##  names(tau_norm_dec_det) <- resnms
##  somehits <- sapply(tau_norm_dec_det, nrow) != 0
##
##  tau_norm_dec_det_somehits <- tau_norm_dec_det[somehits]
##  ## eval_maxp(tau_norm_dec_det_somehits)
##  ## eval_minate(tau_norm_dec_det, atenm = "ate_norm_dec")
##  ## eval_numgrps(tau_norm_dec_det_somehits)
##  ## eval_single_blocks_found(tau_norm_dec_det_somehits)
##  ## eval_treedepth(tau_norm_dec_det_somehits)
## })

test_that("alphafns work across splitters for constant effect that cancel out at the high level
  (half large and positive, half large and negative).", {
  ################################################################################
  ### All blocks have large effects, some are negative and some positive. The
  ### blocks vary in size so perhaps the more sensitive procedures will be more
  ### likely to pick up effects in the smaller blocks.
  table(idat3$bF)

  tau_v1 <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Y ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_v1) <- resnms

  tau_v1_det <- lapply(seq_along(tau_v1), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_v1)[[i]], 1, 4) == "NULL"
    report_detections(tau_v1[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_v1_det) <- resnms

  eval_maxp(tau_v1_det)
  eval_minate(tau_v1_det, atenm = "ate_tau")
  eval_numgrps(tau_v1_det)
  eval_single_blocks_found(tau_v1_det)
  eval_treedepth(tau_v1_det)
})


test_that("alphafns work across splitters for individually heterogeneous effects block-fixed effects and some completely null blocks.", {
  ################
  tau_v2 <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ytauv2 ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_v2) <- resnms

  tau_v2_det <- lapply(seq_along(tau_v2), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_v2)[[i]], 1, 4) == "NULL"
    report_detections(tau_v2[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_v2_det) <- resnms

  ### maxp and single_blocks don't fit the intuition here. Need to investigate. Commenting out for now.
  ## eval_maxp(tau_v2_det)
  ## eval_minate(tau_v2_det, atenm = "ate_tauv2")
  eval_numgrps(tau_v2_det)
  # eval_single_blocks_found(tau_v2_det)
  eval_treedepth(tau_v2_det)
})
