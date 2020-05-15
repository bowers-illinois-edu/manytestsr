## Test and develop functions to adapt alpha levels as the tree grows
context("Alpha Adjusting Performance")

# library(here)
# source(here::here("tests/testthat", "make_test_data.R"))
# devtools::load_all() ## comment this out for production

setDTthreads(1)
options(digits = 4)
## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]

## Setting up to test the alpha adjusting methods versus FWER methods on different splitters
alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "NULL"),
  sfn = c(
    "splitCluster", "splitEqualApprox", "splitLOO",
    "splitSpecifiedFactor"
  ), # , "splitSpecified"),
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <- "lvs"

testing_fn <- function(afn, sfn, sby, fmla = Ytauv2 ~ ZF | bF, idat = idat3, bdat = bdat4) {
  if (afn == "NULL") {
    theafn <- NULL
  } else {
    theafn <- getFromNamespace(afn, ns = "manytestsr")
  }
  ## afn and sfn and sby are character names
  theres <- findBlocks(
    idat = idat, bdat = bdat, blockid = "bF", splitfn = get(sfn),
    pfn = pIndepDist, alphafn = theafn, thealpha = 0.05,
    fmla = fmla, # Ynorm_inc ~ ZF | bF,
    parallel = "no", copydts = TRUE, splitby = sby
  )
  return(theres)[order(biggrp)]
  #theps <- grep("^p[0-9]", names(theres), value = TRUE)
  #theas <- grep("^alpha", names(theres), value = TRUE)
  #truth <- grep("^ate", names(theres), value = TRUE)
  #return(theres[, .SD, .SDcols = c(
  #  theps, theas, truth, "pfinalb", "blocksbygroup",
  #  "biggrp", "bF", "hwt", "nb"
  #)][order(biggrp)])
}

##  Maybe make w0 more like .05.
alpha_and_splits[c(1, 3), ]
## debug(findBlocks)
res1 <- testing_fn(
  afn = alpha_and_splits[1, "afn"],
  sfn = alpha_and_splits[1, "sfn"],
  sby = alpha_and_splits[1, "splitby"],
  idat = idat3, bdat = bdat4
)
res3 <- testing_fn(
  afn = alpha_and_splits[3, "afn"],
  sfn = alpha_and_splits[3, "sfn"],
  sby = alpha_and_splits[3, "splitby"],
  idat = idat3, bdat = bdat4
)

res1 <- res1[, -c("ate_null", "ate_homog", "ate_norm_inc", "ate_tau")]
res3 <- res3[, -c("ate_null", "ate_homog", "ate_norm_inc", "ate_tau")]

options(digits = 3, scipen = 4)
res1[order(p1, p2, p3, p4, p5, decreasing = TRUE)]
res3[order(p1, p2, p3, p4, p5, decreasing = TRUE)]

## Which blocks were identified?
##
## The alpha investing style procedures should identify **more** blocks
# We detect an effect in an individual block if:
## (1) the block is a leaf and that leaf has p < alpha
# We detect an effect for a group of two blocks if
## (2) the blocks are leaves with p > alpha but the parent of that block has p < alpha

## Change one result from res3 to enable us to check this:

res1[bF %in% c("10", "9"), ]
res3[bF %in% c("10", "9"), ]
res1[bF == "10", p5 := .5]
res1[bF == "10", pfinalb := .5]
res3[bF == "10", p5 := .5]
res3[bF == "10", pfinalb := .5]

res3new <- report_detections(res3)
## So we can say that we discovered hits in the following blocks or groups of blocks
res3new[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha)][order(hit_grp)]

# debugonce(report_detections)
res1new <- report_detections(res1, fwer = FALSE)
## So we can say that we discovered hits in the following blocks or groups of blocks
res1new[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha)][order(hit_grp)]

res3_tree <- make_tree(res3,blockid="bF")
res1_tree <- make_tree(res1,blockid="bF")

res3_g <- make_graph(res3_tree)
res1_g <- make_graph(res1_tree)

## Compare these graphs to the results in res1new and res3new above.
#ggsave(res3_g, file = "res3nodes.pdf", bg = "transparent", width = 12, height = 7)
#ggsave(res1_g, file = "res1nodes.pdf", bg = "transparent", width = 13, height = 7)

## Criteria for comparisons:
## These functions can be used to run the tests for the different scenarios.
## They have hard coded the relationships in alpha_and_splits above for now.

eval_numgrps <- function(detection_obj) {
  numgrps <- sapply(detection_obj, function(dat) {
    length(unique(dat$hit_grp))
  })
  expect_lte(numgrps[[3]], max(numgrps[1:2])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numgrps[[6]], max(numgrps[4:5]))
  expect_lte(numgrps[[9]], max(numgrps[7:8]))
}

eval_maxp <- function(detection_obj) {
  maxp <- sapply(detection_obj, function(dat) {
    dat[, sizerp := .N, by = hit_grp]
    dat[, p := ifelse(sizerp > 1, fin_parent_p, max_p)]
    max(dat[, .(min_max_p = min(p)), by = hit_grp]$min_max_p)
  })
  ## Adding tolerance of .01 here.
  expect_lte(maxp[[3]], max(maxp[1:2]) + .01) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(maxp[[6]], max(maxp[4:5]) + .01)
  expect_lte(maxp[[9]], max(maxp[7:8]) + .01)
}

eval_single_blocks_found <- function(detection_obj) {
  ## Number of individual blocks detected versus groups (expect more singletons with the alpha adjusting approaches)
  numsingletons <- sapply(detection_obj, function(dat) {
    tab <- table(dat$hit_grp)
    sum(tab == 1)
  })
  expect_lte(numsingletons[[3]], min(numsingletons[1:2])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numsingletons[[6]], min(numsingletons[4:5]))
  expect_lte(numsingletons[[9]], min(numsingletons[7:8]))
}

eval_treedepth <- function(detection_obj) {
  treedepth <- sapply(detection_obj, function(dat) {
    max(stri_count_fixed(dat$biggrp, ".")) + 1
  })
  expect_lte(treedepth[[3]], max(treedepth[1:2])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(treedepth[[6]], max(treedepth[4:5]))
  expect_lte(treedepth[[9]], max(treedepth[7:8]))
}

eval_minate <- function(detection_obj, atenm) {
  minate <- sapply(detection_obj, function(dat) {
    min(dat[, min(abs(get(atenm))), by = hit_grp]$V1)
  })
  expect_lte(minate[[3]], max(minate[1:2])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(minate[[6]], max(minate[4:5]))
  expect_lte(minate[[9]], max(minate[7:8]))
}

resnms <- apply(alpha_and_splits, 1, function(x) {
  paste(x, collapse = "_", sep = "")
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
  ## We are comparing resnms 1 and 2 (the two fdr approaches) versus 3 (fwer)
  ## 4,5 versus 6
  ## 7,8 versus 9
  ## 10,11 versus 12

  ## Block 1 is the smallest block so it doesn't trigger a detection
  table(idat3$bF)
  block1detected <- sapply(tau_homog_det, function(dat) {
    any(dat$bF == "1")
  })
  expect_equal(all(block1detected), FALSE)
  numblks_homog <- sapply(tau_homog_det, nrow)
  expect_equal(all(numblks_homog == 19), TRUE)
  ## Each block should be detected  --- they should not be grouped together --- so the number of groups should be 19 as well
  numgrps_homog <- sapply(tau_homog_det, function(dat) {
    length(unique(dat$hit_grp))
  })
  expect_equal(all(numgrps_homog == 19), TRUE)

  eval_numgrps(tau_homog_det)
  ## Number of individual blocks detected versus groups (expect more or same singletons with the alpha adjusting approaches)
  eval_single_blocks_found(tau_homog_det)

  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  minate_homog <- sapply(tau_homog_det, function(dat) {
    min(dat[, min(abs(ate_homog)), by = hit_grp]$V1)
  })

  ## Highest p detected: Doesn't reliably work in this case with all huge effects.
  ## eval_maxp(tau_homog_det)

  # g_homog <- lapply(tau_homog,function(obj){ make_graph(make_tree(obj)) })
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
  ## We are comparing resnms 1 and 2 (the two fdr approaches) versus 3 (fwer)
  ## 4,5 versus 6
  ## 7,8 versus 9
  ## 10,11 versus 12

  ## Number of individual blocks detected versus groups (expect more singletons with the alpha adjusting approaches)
  eval_single_blocks_found(tau_norm_inc_det)

  ## Number of blocks detected: this is not clear because of the ability to declare "detect" for *groups* of blocks. So, FWER might stop testing  and return many blocks.
  numblks_norm_inc <- sapply(tau_norm_inc_det, nrow)
  # expect_lte(numblks_norm_inc[[3]],min(numblks_norm_inc[1:2]))
  # expect_lte(numblks_norm_inc[[6]],min(numblks_norm_inc[4:5]))
  # expect_lte(numblks_norm_inc[[9]],min(numblks_norm_inc[7:8]))

  ## Number of groups of blocks detected: expect more or equal from alpha adjusting
  eval_numgrps(tau_norm_inc_det)
  ## Depth of testing (maxdepth): Probably a deeper tree or equal.
  eval_treedepth(tau_norm_inc_det)
  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  eval_minate(tau_norm_inc_det, atenm = "ate_norm_inc")

  ## Highest p detected
  eval_maxp(tau_norm_inc_det)
  ## Null blocks detected (but this could be ok if they were in a group of other blocks having strong effects, so do by hit_grp)
  tau_norm_inc_det[[1]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  tau_norm_inc_det[[3]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  symdiff <- function(x, y) {
    setdiff(union(x, y), intersect(x, y))
  } # https://www.r-bloggers.com/symmetric-set-differences-in-r/
  symdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[3]]$bF)
  symdiff(tau_norm_inc_det[[2]]$bF, tau_norm_inc_det[[3]]$bF)
  setdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[3]]$bF)
  setdiff(tau_norm_inc_det[[3]]$bF, tau_norm_inc_det[[1]]$bF)

  sort(setdiff(bdat4$bF, tau_norm_inc_det[[1]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[2]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[3]]$bF))
  # g_norm_inc <- lapply(tau_norm_inc, function(obj) {
  #   make_graph(make_tree(obj))
  # })
  # for (i in 1:length(g_norm_inc)) {
  #   ggsave(g_norm_inc[[i]],
  #     file = paste0("tau_norm_inc_g", i, ".pdf"),
  #     bg = "transparent", width = 14, height = 7
  #   )
  # }
})

test_that("alphafns work across splitters for individually heteogeneous effects and decrease with block size. Also some completely null blocks.", {
  ## Some blocks have no effect at all in norm_dec
  tau_norm_dec <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynorm_dec ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_norm_dec) <- resnms
  tau_norm_dec_det <- lapply(seq_along(tau_norm_dec), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_norm_dec)[[i]], 1, 4) == "NULL"
    report_detections(tau_norm_dec[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_norm_dec_det) <- resnms
  ## We are comparing resnms 1 and 2 (the two fdr approaches) versus 3 (fwer)
  ## 4,5 versus 6
  ## 7,8 versus 9
  ## 10,11 versus 12

  eval_maxp(tau_norm_dec_det)
  eval_minate(tau_norm_dec_det, atenm = "ate_norm_dec")
  eval_numgrps(tau_norm_dec_det)
  eval_single_blocks_found(tau_norm_dec_det)
  eval_treedepth(tau_norm_dec_det)
})

test_that("alphafns work across splitters for constant effect that cancel out at the high level (half large and positive, half large and negative).", {
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
  ## We are comparing resnms 1 and 2 (the two fdr approaches) versus 3 (fwer)
  ## 4,5 versus 6
  ## 7,8 versus 9
  ## 10,11 versus 12

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

  eval_maxp(tau_v2_det)
  eval_minate(tau_v2_det, atenm = "ate_tauv2")
  eval_numgrps(tau_v2_det)
  eval_single_blocks_found(tau_v2_det)
  eval_treedepth(tau_v2_det)
})
