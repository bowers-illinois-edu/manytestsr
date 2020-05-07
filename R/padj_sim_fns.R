## Functions to assess the procedures

##' Create treatment effect, then add tau and test using the SIUP method
##'
##' This function returns a function that carries with it the environment containing the arguments.
##' The idea is to make parallelization easier.
##' @param idat Individual level data
##' @param trtid Is the name of the treatment numeric, (0,1), variable
##' @param fmla is a simple formula with no blocks --- we are calculating tests within all blocks.
##' @param covariate is the name of a covariate to be used in created covariate dependent treatment effects
##' @return A pvalue for each block
##' @export
padj_test_fn <- function(idat, bdat, blockid, trtid = "trt", fmla = Y ~ trtF | blockF, ybase,
                         prop_blocks_0, tau_fn, tau_size, pfn, afn, p_adj_method = "fdr", nsims, ncores = 1,
                         reveal_and_test_fn, splitfn = NULL, covariate = NULL, splitby = NULL, thealpha=.05) {

    if (!is.null(afn) & is.character(afn)){
        if(afn == "NULL") {
            afn <- NULL
        } else {
            afn <- get(afn)
        }
    }

  datnew <- copy(idat)
  bdatnew <- copy(bdat)

  setkeyv(datnew, blockid)
  setkeyv(bdatnew, blockid)

  set.seed(12345) ## make data the same way for each design. This is a hack.
  datnew$y1new <- create_effects(
    idat = datnew, ybase = ybase, blockid = blockid,
    tau_fn = tau_fn, tau_size = tau_size,
    prop_blocks_0 = prop_blocks_0, covariate = covariate
  )
  datnew[, trueblocks := ifelse((y1new - get(ybase)) == 0, 0, 1)]
  bdat_effects <- datnew[, .(trueblocks = mean(trueblocks)), by = blockid]
  setkeyv(bdat_effects, blockid)
  stopifnot(all.equal(bdatnew[[blockid]], bdat_effects[[blockid]]))
  bdatnew$trueblocks <- bdat_effects$trueblocks

  if (ncores > 1) {
    #p_sims <- future_replicate(nsims, reveal_and_test_fn(
    #  idat = datnew, bdat = bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new",
    #  fmla = fmla, ybase = ybase, prop_blocks_0 = prop_blocks_0,
    #  tau_fn = tau_fn, tau_size = tau_size, pfn = pfn, afn = afn, p_adj_method = p_adj_method, splitfn = splitfn, splitby = splitby, thealpha=thealpha
    #))
  } else {
    p_sims <- replicate(nsims, reveal_and_test_fn(
      idat = datnew, bdat = bdatnew, blockid = blockid, trtid = trtid, y1var = "y1new",
      fmla = fmla, ybase = ybase, prop_blocks_0 = prop_blocks_0,
      tau_fn = tau_fn, tau_size = tau_size, pfn = pfn, afn = afn, p_adj_method = p_adj_method, splitfn = splitfn, splitby = splitby, thealpha=thealpha
    ))
    # ,mc.cores=ncores)
  }
  return(p_sims)
}

##' Repeat experiment, reveal treatment effects from the potential outcomes, test within each block, summarize
##'
##' The function does hypothesis tests within blocks and then summarizes the results  of this testing across the blocks. It is designed for use by padj_test_fn.
##' @param y1var Is the name of the potential outcome to treatment
##' @param p_adj_method Is the argument to p.adjust()
##' @param splitfn Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
##' @param splitby Must be null. Only exists so that we can use the same efffect creation and testing function for both approaches.
##' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
##' @export
reveal_po_and_test <- function(idat, bdat, blockid, trtid, fmla = NULL, ybase, y1var,
                               prop_blocks_0, tau_fn, tau_size, pfn, p_adj_method = "fdr", splitfn = NULL, splitby = NULL, thealpha=.05) {
  stopifnot(is.null(splitfn))
  idat[, newZ := sample(get(trtid)), by = blockid] ## make no effects within block by shuffling treatment
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)] ## reveal relevant potential outcomes.
  idat[, newZF := factor(newZ)] ## the pvalue functions want a factor
  fmla <- Y ~ newZF ## redundant, shouldn't be in arguments if overriding here.
  theps <- idat[, .(p = pfn(fmla = fmla, dat = .SD, parallel = "no")), by = blockid]
  setkeyv(theps, blockid)
  setkeyv(bdat, blockid)
  res <- merge(theps, bdat)
  res[, p_adj := p.adjust(p, method = p_adj_method)]
  fper <- as.numeric(any(res$p <= thealpha))
  fdp <- sum(res$p[res$trueblocks == 0] <= thealpha) / max(1, sum(res$p <= thealpha))
  padj_fdp <- sum(res$p_adj[res$trueblocks == 0] <= thealpha) / max(1, sum(res$p_adj <= thealpha))
  padj_pwr <- sum(res$p_adj[res$trueblocks == 1] <= thealpha) / max(1, sum(res$trueblocks == 1))
  pwr <- sum(res$p[res$trueblocks == 1] <= thealpha) / max(1, sum(res$trueblocks == 1))
  return(c(fper = fper, fdp = fdp, padj_fdp = padj_fdp, padj_pwr = padj_pwr, pwr = pwr))
}

##' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
##'
##' Repeat experiment, reveal treatment effects from the potential outcomes, test within partitions, summarize
##'
##' The function does hypothesis tests within partitions of blocks (including individual blocks depending on the splitting algorithmn) and then summarizes the results  of this testing across the blocks. It very much depends  on padj_test_fn.
##' @param y1var Is the name of the potential outcome to treatment
##' @param p_adj_method Is the argument to p.adjust()
##' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
##' @export
reveal_po_and_test_siup <- function(idat, bdat, blockid, trtid, fmla = Y ~ newZF | blockF, ybase, y1var,
                                    prop_blocks_0, tau_fn, tau_size, pfn, afn, p_adj_method = NULL, copydts = FALSE, splitfn, splitby, thealpha=.05) {
    if (!is.null(afn) & is.character(afn)){
        if(afn == "NULL") {
            afn <- NULL
        } else {
            afn <- get(afn)
        }
    }
  ## 	message(tau_size)
  ## 	make relationship between treatment and outcome 0 on average in preparation for power analysis and error rate sims
  idat[, newZ := sample(get(trtid)), by = blockid]
  ## Then reveal an observed outcome that contains a treatment effect from y1var as a function of newZ. So the treatment effect is known.
  idat[, Y := get(y1var) * newZ + get(ybase) * (1 - newZ)]
  ## bdat in findBlocks doesn't really  carry important information other than block id
  setkeyv(bdat, blockid)
  idat[, newZF := factor(newZ)]
  fmla <- as.formula(paste("Y~newZF|", blockid, sep = ""))
  res <- findBlocks(
    idat = idat, bdat = bdat, pfn = pfn, alphafn=afn, splitfn = splitfn, thealpha = thealpha,
    fmla = fmla, parallel = "no", blockid = blockid, copydts = copydts, splitby = splitby
  )
  setkeyv(res, blockid)
  res_report<- report_detections(res, fwer = is.null(afn), only_hits = FALSE)
  if(all(is.na(res_report$hit_grp))){
    return( as.matrix(t(c(fper=0,fdp=0,tdp=0,ngrps=1))) )
    } else {
  ##  Now aggregate to the partition  (or test) rather than the block
  res2 <- res[, .(p = unique(pfinalb), truetau = as.numeric(any(trueblocks == 1))), by = biggrp]
  ## Summarize errors across the splits / blocks
  res3 <- res2[, .(
    fper = as.numeric(any(p <= thealpha)),
    fdp = sum(as.numeric(p <= thealpha) * (1 - truetau)) / max(1, sum(as.numeric(p <= thealpha))),
    tdp = sum(as.numeric(p <= thealpha) * truetau) / max(1, sum(truetau)),
    ngrps = .N
  )]

  return(as.matrix(res3))
    }
}


##' Calculate the error and success proportions of tests for a single iteration
##'
##' This function takes output from [findBlocks] where a or equivalent bottom-up testing function [bottom_up]
##' and returns the proportions of errors made. This means that the input to findBlocks includes a column containing a true block-level effect.
##' Repeated uses of this function allow us to assess false discovery rates and family wise error rates among other metrics of testing success.
##' @param testobj Is an object arising from [findBlocks] or [bottom_up].
##' @param p_adj_method Is the argument to p.adjust()
##' @return False positive proportion out of the tests across the blocks, The false discovery rate (proportion rejected of false nulls out of all rejections), the power of the adjusted tests across blocks (the proportion of correctly rejected hypotheses out of all correct hypotheses --- in this case correct means non-null), and power of the unadjusted test (proportion correctly rejected out of  all correct hypothesis, but using unadjusted p-values).
##' @export
calc_errs <- function(testobj,
                      truevar_name,
                      trueeffect_tol=.Machine$double.eps,
                      thealpha=.05
                      ){

    testobj <- copy(testobj)

    if(length(grep("biggrp",names(testobj)))>0){
        ## this is for the top-down/split and test method
        detobj <- report_detections(testobj, fwer=FALSE, only_hits=FALSE)
        detnodes  <- detobj[, .(hit=as.numeric(unique(hit)),
                                anynull = as.numeric(any( abs(get(truevar_name)) <= trueeffect_tol)),
                                allnull = as.numeric(all( abs(get(truevar_name)) <= trueeffect_tol)),
                                anynotnull = as.numeric(any( abs(get(truevar_name)) > trueeffect_tol ))
                                ), by = hit_grp]

    } else {
        ## This is for the bottom-up/test every block method
        testobj[,hit:=max_p < thealpha]

        detnodes <- testobj[, .(hit = as.numeric(hit),
                                anynull = as.numeric(abs(get(truevar_name)) <= trueeffect_tol),
                                allnull = as.numeric(abs(get(truevar_name)) <= trueeffect_tol),
                                anynotnull = as.numeric(abs(get(truevar_name)) > trueeffect_tol)
                                )]
    }

    deterrs <- detnodes[,.(nreject=sum(hit),
                           naccept=sum(1-hit),
                           prop_reject = mean(hit),
                           prop_accept = mean(1-hit), # or 1-prop_reject
                           ## rejections of false null /detections of true non-null
                           ## (if any of the blocks have an effect, then we count this as a correct rejection or correct detection)
                           true_pos_prop = mean(hit*anynotnull),
                           ## If we reject but *all* of the blocks have no effects, this is a false positive error
                           ## If we reject but only one of them has no effects but the other has effects,
                           ## then this is not an error --- but a correct detection
                           false_pos_prop = mean(hit*allnull),
                           ## If we do not reject and all blocks are truly null, then we have no error.
                           true_neg_prop = mean((1-hit)*allnull),
                           ## If we do not reject/detect and at least one of the blocks actually has an effect, we have
                           ## a false negative error --- a failure to detect the truth
                           false_neg_prop = mean((1-hit)*anynotnull),
                           ## Now look at false and true discoveries: false rejections as a proportion of rejections
                           false_disc_prop = sum(hit*allnull)/max(1,sum(hit)),
                           true_disc_prop = sum(hit*anynotnull)/max(1,sum(hit)),
                           true_nondisc_prop = sum((1-hit)*allnull)/max(1,sum(1-hit)),
                           false_nondisc_prop = sum((1-hit)*anynotnull)/max(1,sum(1-hit))
                           )]

    return(deterrs)
}

