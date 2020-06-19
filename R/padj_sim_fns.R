# Functions to assess the procedures via simulation. Mostly for use in the paper itself.

#' Create treatment effect, then add tau and test using the SIUP method
#'
#' This function returns a function that carries with it the environment containing the arguments.
#' The idea is to make parallelization easier.
#' @param idat Individual level data
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  | block where treatment assignment and block must be factors.
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param nsims Is the number of simulations to run --- each simulation uses the same treatment effects be re-assigns treatment (re-shuffles treatment and re-reveals the observed outcomes as a function of the potential outcomes)
#' @param ncores Tells p-value functions how many cores to use. Mostly ignored in use of this function because we are tending to parallelize at higher loops.
#' @param splitfn A function to split the data into two pieces --- using bdat
#' @param covariate is the name of a covariate to be used in created covariate dependent treatment effects
#' @param splitby A string indicating which column in bdat contains a variable to guide splitting (for example, a column with block sizes or block harmonic mean weights or a column with a covariate (or a function of covariates))
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @return A pvalue for each block
#' @export
padj_test_fn <- function(idat, bdat, blockid, trtid = "trt", fmla = Y ~ trtF | blockF, ybase,
                         prop_blocks_0, tau_fn, tau_size, pfn, afn, p_adj_method, nsims, ncores = 1,
                         splitfn = NULL, covariate = NULL, splitby = NULL, thealpha = .05) {
  if (!is.null(afn) & is.character(afn)) {
    if (afn == "NULL") {
      afn <- NULL
    } else {
      afn <- get(afn)
    }
  }

  datnew <- copy(idat)
  bdatnew <- copy(bdat)

  setkeyv(datnew, blockid)
  setkeyv(bdatnew, blockid)

  set.seed(12345) # make data the same way for each design. This is a hack to have the seed within the function.
  datnew$y1new <- create_effects(
    idat = datnew, ybase = ybase, blockid = blockid,
    tau_fn = tau_fn, tau_size = tau_size,
    prop_blocks_0 = prop_blocks_0, covariate = covariate
  )
  datnew[, trueblocks := ifelse(abs(y1new - get(ybase)) <= .Machine$double.eps, 0, 1)]
  datnew[, trueate := mean(y1new - get(ybase)), by = blockid]
  bdat_effects <- datnew[, .(
    trueblocks = unique(trueblocks),
    trueate = mean(trueate)
  ), by = blockid]
  stopifnot(all(bdat_effects$trueblocks %in% c(0, 1)))
  setkeyv(bdat_effects, blockid)
  bdatnew <- bdatnew[bdat_effects, ]
  # bdatnew$trueblocks <- bdat_effects$trueblocks

  if (p_adj_method == "split") {
    reveal_and_test_fn <- reveal_po_and_test_siup
  } else {
    reveal_and_test_fn <- reveal_po_and_test
  }

  if (ncores > 1) {
    message("Not parallelizing this part of the loop.")
    # 	 p_sims <- future_replicate(nsims, reveal_and_test_fn(
    # 	  idat = datnew, bdat = bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new",
    # 	  fmla = fmla, ybase = ybase, prop_blocks_0 = prop_blocks_0,
    # 	  tau_fn = tau_fn, tau_size = tau_size, pfn = pfn, afn = afn, p_adj_method = p_adj_method, splitfn = splitfn, splitby = splitby, thealpha=thealpha
    # 	 ))
  } else {
    p_sims_lst <- replicate(nsims, reveal_and_test_fn(
      idat = datnew, bdat = bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new",
      fmla = fmla, ybase = ybase, prop_blocks_0 = prop_blocks_0,
      tau_fn = tau_fn, tau_size = tau_size, pfn = pfn, afn = afn, p_adj_method = p_adj_method, splitfn = splitfn, splitby = splitby, thealpha = thealpha
    ), simplify = FALSE)
    p_sims <- data.table::rbindlist(p_sims_lst)
  }
  return(p_sims)
}

#' Repeat experiment, reveal treatment effects from the potential outcomes, test within each block, summarize
#'
#' The function does hypothesis tests within blocks and then summarizes the results  of this testing across the blocks. It is designed for use by padj_test_fn.
#' @param idat Individual level data
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  where treatment assignment is a factor
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param p_adj_method Is an input to p.adjust like "fdr" or "holm" to use adjust_block_tests
#' @param afn NULL. Not used. Included here to enable same functions to call this function and reveal_po_and_test_siup
#' @param splitfn Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
#' @param splitby Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @param copydts TRUE or FALSE. TRUE if using findBlocks standalone. FALSE if copied objects are being sent to findBlocks from other functions.
#' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
#' @export
reveal_po_and_test <- function(idat, bdat, blockid, trtid, fmla = NULL, ybase, y1var,
                               prop_blocks_0, tau_fn, tau_size, pfn, p_adj_method = "fdr",
                               afn = NULL, splitfn = NULL, splitby = NULL, thealpha = .05, copydts = FALSE) {
  stopifnot(is.null(splitfn) | splitfn == "NULL")
  idat[, newZ := sample(get(trtid)), by = blockid] # make no effects within block by shuffling treatment, this is the engine of variability in the sim
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)] # reveal relevant potential outcomes with possible known effect
  idat[, newZF := factor(newZ)] # the pvalue functions want a factor
  fmla <- Y ~ newZF
  # idat[, Y := get(y1var) * get(trtid) + get(ybase) * (1 - get(trtid))] # reveal relevant potential outcomes with possible known effect

  res <- adjust_block_tests(
    idat = idat, bdat = bdat, blockid = blockid, p_adj_method = p_adj_method,
    pfn = pfn, fmla = fmla,
    copydts = copydts
  )
  res[, hit := max_p < thealpha] # doing FDR==.05 for now. All that really matters is some fixed number.

  errs <- calc_errs(
    testobj = res,
    truevar_name = "trueate",
    trueeffect_tol = .Machine$double.eps,
    blockid = blockid
  )

  return(errs)
}

#' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
#'
#' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
#'
#' The function does hypothesis tests within partitions of blocks (including individual blocks depending on the splitting algorithmn) and then summarizes the results  of this testing across the blocks. It very much depends  on padj_test_fn.
#' @param idat Data at the unit level.
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param trtid Is the name of the treatment numeric, (0,1), variable
#' @param fmla A formula with outcome~treatment assignment  | block where treatment assignment and block must be factors.
#' @param ybase Is the potential outcome to control upon which the treatment effect will be built
#' @param y1var Is the name of the potential outcome to treatment
#' @param prop_blocks_0 Is the proportion of blocks with no effects at all
#' @param tau_fn Is a function that turns ybase into the potential outcome under treatment --- it is a treatment effect creating function.
#' @param tau_size Is the parameter for the tau_fn --- like the true average effect size within a block.
#' @param pfn A function to produce pvalues --- using idat.
#' @param afn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
#' @param p_adj_method Must be "split" here.
#' @param copydts TRUE or FALSE. TRUE if using findBlocks standalone. FALSE if copied objects are being sent to findBlocks from other functions.
#' @param splitfn A function to split the data into two pieces --- using bdat
#' @param splitby A string indicating which column in bdat contains a variable to guide splitting (for example, a column with block sizes or block harmonic mean weights or a column with a covariate (or a function of covariates))
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
#' @export
reveal_po_and_test_siup <- function(idat, bdat, blockid, trtid, fmla = Y ~ newZF | blockF, ybase, y1var,
                                    prop_blocks_0, tau_fn, tau_size, pfn, afn, p_adj_method = "split",
                                    copydts = FALSE, splitfn, splitby, thealpha = .05) {
  if (!is.null(afn) & is.character(afn)) {
    if (afn == "NULL") {
      afn <- NULL
    } else {
      afn <- get(afn)
    }
  }
  # 	message(tau_size)
  # 	make relationship between treatment and outcome 0 on average in preparation for power analysis and error rate sims
  idat[, newZ := sample(get(trtid)), by = blockid]
  # Then reveal an observed outcome that contains a treatment effect from y1var as a function of newZ. So the treatment effect is known.
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]
  # bdat in findBlocks doesn't really  carry important information other than block id
  # setkeyv(bdat, blockid)
  idat[, newZF := factor(newZ)]
  fmla <- as.formula(paste("Y~newZF|", blockid, sep = ""))
  # idat[, Y := get(y1var) * get(trtid) + get(ybase) * (1 - get(trtid))] # reveal relevant potential outcomes with possible known effect

  res <- findBlocks(
    idat = idat, bdat = bdat, blockid = blockid, splitfn = splitfn,
    pfn = pfn, alphafn = afn, thealpha = thealpha,
    fmla = fmla,
    parallel = "no", copydts = copydts, splitby = splitby
  )

  errs <- calc_errs(
    testobj = res,
    truevar_name = "trueate",
    trueeffect_tol = .Machine$double.eps,
    blockid = blockid
  )

  return(errs)
}


#' Calculate the error and success proportions of tests for a single iteration
#'
#' This function takes output from [findBlocks] where a or equivalent bottom-up testing function [bottom_up]
#' and returns the proportions of errors made. This means that the input to findBlocks includes a column containing a true block-level effect.
#' Repeated uses of this function allow us to assess false discovery rates and family wise error rates among other metrics of testing success.
#' @param testobj Is an object arising from [findBlocks] or [bottom_up].
#' @param truevar_name Is a string indicating the name of the variable containing the true underlying causal effect (at the block level).
#' @param trueeffect_tol Is the smallest effect size below which we consider the effect to be zero (by default is it floating point zero).
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)
#' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
#' @export
calc_errs <- function(testobj,
                      truevar_name,
                      trueeffect_tol = .Machine$double.eps,
                      blockid = "bF",
                      thealpha = .05) {
  simp_summary <- function(x) {
    list(min(x, na.rm = TRUE), mean(x, na.rm = TRUE), median(x, na.rm = TRUE), max(x, na.rm = TRUE))
  }

  if (length(grep("biggrp", names(testobj))) > 0) {
    # this is for the top-down/split and test method
    detobj <- report_detections(testobj, fwer = FALSE, only_hits = FALSE)
    detobj[, hit := as.numeric(detobj$hit)]

    detnodes <- detobj[, .(
      hit = unique(hit),
      anynull = as.numeric(any(abs(get(truevar_name)) <= trueeffect_tol)),
      allnull = as.numeric(all(abs(get(truevar_name)) <= trueeffect_tol)),
      anynotnull = as.numeric(any(abs(get(truevar_name)) > trueeffect_tol)),
      grpsize = .N
    ), by = hit_grp]
    setkey(detnodes, hit_grp)
    # prop_hits <- mean(detobj$hit)

    detnodes_effects <- detobj[, simp_summary(get(truevar_name)), by = list(hit, hit_grp)]
    setnames(detnodes_effects, c("hit", "hit_grp", "minate", "meanate", "medianate", "maxate"))
    # detnodes_effects <- cbind(detnodes_effects1[hit==1,],detnodes_effects1[hit==0,.(min0=minate,mean0=meanate,median0=medianate,max0=maxate)])
    setkey(detnodes_effects, hit_grp)
    detnodes <- detnodes[detnodes_effects, ]
  } else {
    # This is for the bottom-up/test every block method
    detobj <- testobj[, .(
      blockid = get(blockid),
      hit = max_p < thealpha,
      hit_grp = get(blockid),
      truevar_name = get(truevar_name)
    )]
    setnames(detobj, "truevar_name", truevar_name)

    detnodes <- detobj[, .(
      hit = as.numeric(hit),
      anynull = as.numeric(abs(get(truevar_name)) <= trueeffect_tol),
      allnull = as.numeric(abs(get(truevar_name)) <= trueeffect_tol),
      anynotnull = as.numeric(abs(get(truevar_name)) > trueeffect_tol),
      trueeffect = get(truevar_name),
      grpsize = 1,
      blockid = blockid
    )]
    setkey(detnodes, blockid)
    detnodes_effects <- detobj[, as.list(simp_summary(get(truevar_name))), by = list(hit, hit_grp)]
    setnames(detnodes_effects, c("hit", "hit_grp", "minate", "meanate", "medianate", "maxate"))
    # detnodes_effects <- cbind(detnodes_effects1[hit==1,],detnodes_effects1[hit==0,.(min0=minate,mean0=meanate,median0=medianate,max0=maxate)])
    setkey(detnodes_effects, hit_grp)
    detnodes <- detnodes[detnodes_effects, ]
  }
  detates1 <- detnodes[, .(
    minate = min(minate),
    meanate = mean(meanate),
    medate = median(medianate),
    maxate = max(maxate)
  ), by = hit]
  if (nrow(detates1) == 1) {
    detates2 <- rbind(detates1, list(1 - detates1$hit, NA, NA, NA, NA))
  }

  detates <- do.call("cbind", lapply(c(1, 0), function(i) {
    obj <- detates1[hit == i, ]
    setnames(obj, paste0(names(obj), i))
    return(obj)
  }))

  deterrs <- detnodes[, .(
    nreject = sum(hit),
    naccept = sum(1 - hit),
    prop_reject = mean(hit),
    prop_accept = mean(1 - hit), # or 1-prop_reject
    # rejections of false null /detections of true non-null
    # (if any of the blocks have an effect, then we count this as a correct rejection or correct detection)
    true_pos_prop = mean(hit * anynotnull),
    # If we reject but *all* of the blocks have no effects, this is a false positive error
    # If we reject but only one of them has no effects but the other has effects,
    # then this is not an error --- but a correct detection
    false_pos_prop = mean(hit * allnull),
    # If we do not reject and all blocks are truly null, then we have no error.
    true_neg_prop = mean((1 - hit) * allnull),
    # If we do not reject/detect and at least one of the blocks actually has an effect, we have
    # a false negative error --- a failure to detect the truth
    false_neg_prop = mean((1 - hit) * anynotnull),
    # Now look at false and true discoveries: false rejections as a proportion of rejections
    false_disc_prop = sum(hit * allnull) / max(1, sum(hit)),
    true_disc_prop = sum(hit * anynotnull) / max(1, sum(hit)),
    true_nondisc_prop = sum((1 - hit) * allnull) / max(1, sum(1 - hit)),
    false_nondisc_prop = sum((1 - hit) * anynotnull) / max(1, sum(1 - hit)),
    meangrpsize = mean(grpsize),
    medgrpsize = median(grpsize)
  )]
  # One row of results at the end.
  detresults <- cbind(deterrs, detates)
  return(detresults)
}
