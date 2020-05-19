# Functions for combining splitting and testing

##' Test, Split, Repeat
##'
##' Split and test.
##' @param idat Data at the unit level.
##' @param bdat Data at the block level.
##' @param splitfn A function to split the data into two pieces --- using bdat
##' @param pfn A function to produce pvalues --- using idat.
##' @param alphafn A function to adjust alpha at each step. Takes one or more p-values plus a stratum or batch indicator.
##' @param sims Number of permutations for permutation-based testing
##' @param copydts TRUE or FALSE. TRUE if using findBlocks standalone. FALSE if copied objects are being sent to findBlocks from other functions.
##' @param splitby A string indicating which column in bdat contains a variable to guide splitting (for example, a column with block sizes or block harmonic mean weights or a column with a covariate (or a function of covariates))
##' @param blocksize A string with the name of the column in bdat contains information about the size of the block (or other determinant of the power of tests within that block, such as harmonic mean weight of the block or variance of the outcome within the block.)
##' @return A data.table containing information about the sequence of splitting and testing
##' @importFrom stringi stri_count_fixed stri_split_fixed stri_split stri_sub stri_replace_all stri_extract_last
##' @importFrom digest digest getVDigest
##' @export
findBlocks <- function(idat, bdat, blockid = "block", splitfn, pfn, alphafn = NULL, simthresh = 20,
                       sims = 1000, maxtest = 100, thealpha = 0.05,
                       fmla = YContNorm ~ trtF | blockF,
                       parallel = "multicore", copydts = FALSE, splitby = "hwt", blocksize="hwt") {
  if (copydts) {
    bdat <- data.table::copy(bdat)
    idat <- data.table::copy(idat)
  }

  ## Step 1: Root of tree. One test. We are numbering simply following literature on complete or perfect binary trees
  i <- 1L
  bdat[, p1 := pfn(
    fmla = fmla, dat = idat,
    simthresh = simthresh, sims = sims, parallel = parallel
  )]
  bdat[, pfinalb := p1] ## just to initialize pfinalb
  ## Record the names of the group and pvalues.
  gnm <- "g1"
  pnm <- "p1"
  bdat[, biggrp := gl(1, nrow(bdat), labels = "1")] ## initialize the biggrp factor var
  bdat[, g1 := biggrp]
  if (is.null(alphafn)) {
    bdat[, alpha1 := thealpha]
  } else {
    ## Current thought: the FDR controlling procedure should not be more conservative
    ## than the FWER controlling procedure.
    ## So, use which ever alpha is bigger at the overall root.
    bdat[, alpha1 := max(alphafn(pval = unique(p1), batch = unique(g1), nodesize=sum(get(blocksize))), thealpha)]
  }
  ## pbprev <- bdat[,.(p=unique(p1),biggrp=unique(g1),alpha1=unique(alpha1),batch="p1")]
  # If p1 > alpha1 then stop testing and return the data.
  bdat[, testable := (pfinalb < alpha1)]
  pbtracker <- data.table(
    p = unique(bdat$p1), biggrp = factor("1"),
    alpha1 = unique(bdat$alpha1), batch = "p1",
    testable = unique(bdat$testable),
    nodenum = "1",
    depth = 1L,
    nodesize=sum(bdat[,get(blocksize)])
  )
  bdat[, nodenum_current := "1"] 
  bdat[, nodenum_prev:= "0"] 
  ## there is only one test at this root level.
  ## stopifnot(length(unique(dat$testable))==1)
  setkeyv(idat, blockid)
  setkeyv(bdat, "testable")
  bdat[, blocksbygroup := length(unique(get(blockid))), by = biggrp]
  ## samep <- FALSE ## initialize tracking of p-values to check if they begin to be repeated
  if (!all(bdat$testable)) {
    return(bdat)
  }
  nodeidfn <- function(d){
      crc32hash <- getVDigest(algo = "crc32")
      crc32hash(d)
  }
  # Step 2: Iterate through the tree
  ## Keep testing until no testing possible: criteria all blocks all size 1 and or
  ## all p_i > alpha_i or there is no change in the final p-values OR simulation
  ## limits reached (for testing of the algorithm).

  while (any(bdat$testable, na.rm = TRUE) & i < maxtest) {
   ##   if(i==10){ browser() }
    i <- i + 1L
    # if (i == 1 | i > 10) {  message("Split number: ", i) }
    gnm <- paste0("g", i) ## name of the grouping variable for the current split
    pnm <- paste0("p", i) ## name of the p-value variable for the current split
    alphanm <- paste0("alpha", i)
    ## Set Group to NA for blocks where we have to stop testing
    bdat[(testable), (gnm) := splitfn(bid = get(blockid), x = get(splitby)), by = biggrp] ## get(paste0('g',i-1))]
    ## maybe improve this next with
    ## https://stackoverflow.com/questions/33689098/interactions-between-factors-in-data-table
    if (i == 2) {
      bdat[(testable), nodenum_prev := nodenum_current]
      #bdat[(testable), nodenum_current := fifelse(get(gnm) == "0", nodenum_prev * 2L , ( nodenum_prev * 2L )  + 1L)]
      bdat[(testable), nodenum_current := fifelse(get(gnm) == "0", nodeidfn( paste0(nodenum_prev,"0") ) ,
                                                  nodeidfn(paste0(nodenum_prev,"1")))]
    } else {
      bdat[(testable), nodenum_prev := nodenum_current]
      ## bdat[(testable), nodenum_current := fifelse(get(gnm) == "0", nodenum_prev * 2L , ( nodenum_prev * 2L )  + 1L), by = biggrp]
      bdat[(testable), nodenum_current := fifelse(get(gnm) == "0", nodeidfn( paste0(nodenum_prev,"0") ) ,
                                                  nodeidfn(paste0(nodenum_prev,"1"))), by = biggrp]
      ## How to test/efficiently or periodically for duplicated nodenums?
    }
    bdat[(testable), biggrp := interaction(biggrp, nodenum_current, drop = TRUE)]
    bdat[, biggrp := droplevels(biggrp)] ## annoying to need this extra step given drop=TRUE
    bdat[, nodesize:=sum(get(blocksize)),by=biggrp]
    ## Now merge idat and bdat again since the test has to be at the idat level
    idat[bdat, c("testable", "biggrp") := mget(c("i.testable", "i.biggrp")), on = blockid]
    idat[, biggrp := droplevels(biggrp)] ## annoying to need to do this
    pb <- idat[(testable), list(p = pfn(
      fmla = fmla, dat = .SD, simthresh = simthresh,
      sims = sims, parallel = parallel
    )), by = biggrp]
    pb[, depth := i]
    ## This next could be made more efficient without string splitting
    pb[, nodenum := stri_split_fixed(biggrp, ".", simplify = TRUE)[, i]]
    ## call "blocksize" the sum of the block sizes within group
    pb[bdat, nodesize:=i.nodesize, on="biggrp"]
    pbtracker <- rbind(pbtracker[, .(p, biggrp, batch, testable, nodenum, depth, nodesize)],
      pb[, batch := pnm],
      fill = TRUE
    )
    setnames(pb, "p", pnm)
    setkeyv(pb, "biggrp")
    bdat[pb, (pnm) := get(paste0("i.", pnm)), on = "biggrp"]
    ## Recall that the siup requires a hierarchy of p-values. So, only
    ## keep the maximum p-value associated with any given block.
    ## bdat[, pfinalbminus1 := pfinalb] ## keep track of previous pvalues
    bdat[(testable), pfinalb := pmax(get(pnm), pfinalb)] # get(paste0("p", i - 1)))]
    ## If the final p-values are not changing then we stop. Maybe change this if alphafn is provided.
    # samep <- identical(bdat$pfinalb, bdat$pfinalbminus1) ## see top of the while() statement
    ## If the p<alpha, continue testing and splitting unless there is only one block
    ## left.
    ## Now decide which blocks (units) can be tested again.
    ## If a split contains only one block. We cannot test further.
    bdat[, blocksbygroup := length(unique(get(blockid))), by = biggrp]
    ## Now update alpha if we are trying to control FDR or mFDR rather than FWER
    if (is.null(alphafn)) {
      bdat[, (alphanm) := thealpha]
      ## Recall that := **updates** values. So, it doesn't overwrite (testable==FALSE) values
      bdat[, testable := fifelse((pfinalb <= get(alphanm)) & (blocksbygroup > 1), TRUE, FALSE)]
    } else {
      if (i == 2) {
        pbtracker[, (alphanm) := alphafn(pval = p, batch = batch, nodesize=nodesize)]
      } else {
        mid_roots <- pbtracker[depth == (i - 1) & (testable), nodenum]
        find_paths <- function(){
          tmp <- pbtracker[depth==i,]
          tmp[,leaves:= stri_extract_last(as.character(biggrp),regex="\\.[:alnum:]*$")]
          tmp[,paths:= stri_replace_all(as.character(biggrp),replacement="",fixed=leaves)]
          path_dat <-tmp[,.(thepath=paste(unique(paths),paste(leaves,sep="",collapse=""),sep="")),by=paths]
          path_vec <- stri_split_fixed(path_dat$thepath,".") #,simplify=TRUE)#[1,]
          ##stopifnot(all(path_vec %in% pbtracker$nodenum))
          return(path_vec)
         }
        setkey(pbtracker, nodenum)
        thepaths <- find_paths()
        for (j in 1:length(thepaths)){
          pbtracker[J(thepaths[[j]]), (alphanm) := alphafn(pval = p, batch = depth, nodesize = nodesize)]
        }
      }
      setkey(pbtracker, biggrp)
      bdat[pbtracker, (alphanm) := get(paste0("i.", alphanm)), on = "biggrp"]
      ## Recall that := **updates** values. So, it doesn't overwrite
      ## (testable==FALSE) values. Here we use current p rather than max of previous p.
      bdat[, testable := fifelse((get(pnm) <= get(alphanm)) & (blocksbygroup > 1), TRUE, FALSE)]
      testable_grps <- bdat[!is.na(testable), .(testable = unique(testable)), keyby = "biggrp"]
      pbtracker[testable_grps, testable := get("i.testable")]
    }
    setkeyv(bdat, "testable") ## for binary search speed
    ## message(paste(unique(signif(bdat$pfinalb,4)),collapse=' '),appendLF = TRUE)
    ## message('Number of blocks left to test: ', sum(bdat$testable))
  }
  return(bdat)
  ## }
}
