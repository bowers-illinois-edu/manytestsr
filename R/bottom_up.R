## Functions for bottom-up testing

##' Test every block and adjust the p-values
##'
##'
##' @param idat Data at the unit level.
##' @param bdat Data at the block level.
##' @param pfn A function to produce pvalues --- using idat.
##' @param sims Number of permutations for permutation-based testing
##' @param copydts TRUE or FALSE. TRUE if using standalone. Could be FALSE if copied objects are being sent to this function from other functions. Basically the question is whether we want this function to change the block or individual level datasets in the main environment or return a new object.
##' @return A data.table containing information about the sequence tests
##' @importFrom stringi stri_count_fixed stri_split_fixed stri_split stri_sub
##' @export
adjust_block_tests <- function(idat, bdat, blockid = "block", pfn, p_adj_method, simthresh = 20,
                               sims = 1000, fmla = YContNorm ~ trtF | blockF,
                               parallel = "multicore", copydts = FALSE) {
    if (copydts) {
        bdat <- data.table::copy(bdat)
        idat <- data.table::copy(idat)
    }
    theps <- idat[, .(p = pfn(fmla = fmla, dat = .SD, parallel = "no")), by = blockid]
    setkeyv(theps, blockid)
    setkeyv(bdat, blockid)
    res <- merge(theps, bdat)
    res[, max_p := p.adjust(p, method = p_adj_method)]
    return(res)
}


