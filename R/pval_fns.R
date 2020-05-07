# Functions for producing p-values

##' P-value function: T-test
##'
##' These functions accept a data frame and perhaps test specific arguments
##' (like whether or not the test will be asympotic or simulation based). It
##' produces a p-value.
##'
##' @param dat An object inheriting from class data.frame
##' @param sims Either NULL (meaning use an asymptotic reference dist) or a
##' number (meaning sampling from the randomization distribution implied by the
##' formula)
##' @return A p-value
##' @importFrom coin oneway_test pvalue approximate exact asymptotic
##' @export
pOneway <- function(dat, fmla = YContNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                    parallel = "no", ncpu = NULL) {
  require(coin)
  theresponse <- all.vars(fmla)[attr(terms(fmla), "response")]
  ## if we get a constant outcome, then p=1.
  ## no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  if (is.null(simthresh) | nrow(dat) > simthresh) {
    thep <- pvalue(oneway_test(fmla, data = dat))[[1]]
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- parallel::detectCores()
    }
    if (parallel == "no") {
      ncpu <- 1
    }
    thep <- pvalue(oneway_test(fmla,
      data = dat,
      distribution = approximate(
        nresample = sims,
        parallel = parallel, ncpus = ncpu
      )
    ))[[1]]
  }
  return(as.numeric(thep))
}

##' P-value function: Wilcox Test
##'
##' These functions accept a data frame and perhaps test specific arguments
##' (like whether or not the test will be asympotic or simulation based). It
##' produces a p-value.
##'
##' @param dat An object inheriting from class data.frame
##' @param sims Either NULL (meaning use an asymptotic reference dist) or a
##' number (meaning sampling from the randomization distribution implied by the
##' formula)
##' @return A p-value
##' @importFrom coin wilcox_test pvalue approximate exact asymptotic
##' @export
pWilcox <- function(dat, fmla = YContNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                    parallel = "no", ncpu = NULL) {
  require(coin)
  theresponse <- all.vars(fmla)[attr(terms(fmla), "response")]
  ## if we get a constant outcome, then p=1.
  ## no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  if (is.null(simthresh) | nrow(dat) > simthresh) {
    thep <- pvalue(wilcox_test(fmla, data = dat))[[1]]
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- parallel::detectCores()
    }
    if (parallel == "no") {
      ncpu <- 1
    }
    thep <- pvalue(wilcox_test(fmla,
      data = dat,
      distribution = approximate(
        nresample = sims,
        parallel = parallel, ncpus = ncpu
      )
    ))[[1]]
  }
  return(as.numeric(thep))
}



##' P-value function: Disco
##'
##' These functions accept a data frame and perhaps test specific arguments
##' (like whether or not the test will be asympotic or simulation based). It
##' produces a p-value.
##'
##' @param dat An object inheriting from class data.frame
##' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment*block
##' @param sims Either NULL (meaning use an asymptotic reference dist) or a
##' number (meaning sampling from the randomization distribution implied by the
##' formula)
##' @param  index is a real number in  (0,2] to transform the inter-unit distances
##' @return A p-value
##' @importFrom coin wilcox_test pvalue approximate exact asymptotic
##' @export
pAnova <- function(dat, fmla = YContNorm ~ trtF * blockF, simthresh = 20, sims = 1000,
                   parallel = "no", ncpu = NULL, index = 1, blocks = blockF) {
  require(energy)
  theresponse <- all.vars(fmla)[attr(terms(fmla), "response")]
  ## if we get a constant outcome, then p=1.
  ## no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  if (is.null(simthresh) | nrow(dat) > simthresh) {
    thep <- anova1$`Pr(>F)`
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- parallel::detectCores()
    }
    if (parallel == "no") {
      ncpu <- 1
    }
    thep <- pvalue(wilcox_test(fmla,
      data = dat,
      distribution = approximate(
        nresample = sims,
        parallel = parallel, ncpus = ncpu
      )
    ))[[1]]
  }
  return(as.numeric(thep))
}


##' P-value function: Independence Treatment Distance Test
##'
##' These functions accept a data frame and perhaps test specific arguments
##' (like whether or not the test will be asympotic or simulation based). It
##' produces a p-value.
##'
##' For now, this  function  does an omnibus-style chi-square  test using (1)
##' the ratio  of  distances to controls to distances to treated observations
##' within block;  (2)  the  rank of  distances to  controls for each unit; and
##' (3) the raw outcome.
##' @param dat An object inheriting from class data.frame
##' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
##' @param sims Either NULL (meaning use an asymptotic reference dist) or a
##' number (meaning sampling from the randomization distribution implied by the
##' formula)
##' @param simthresh is the size of the data below which we use direct permutations for p-values
##' @param  groups is a vector defining the groups within which the inter-unit  distances are calculated. Not used here.
##' @param  distfn is  a function that produces one or more vectors (a data frame or matrix) of the  same number of  rows as the dat
##' @param parallel is "no" then parallelization is not required, otherwise it is "multicore" or "snow" (see  help for coin::approximate()
##' @param ncpu is number of cpus  to be used for parallel operation.
##' @return A p-value
##' @importFrom coin independence_test pvalue approximate exact asymptotic
##' @importFrom Rfast Rank
##' @export
pIndepDist <- function(dat, fmla = YcontNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                       parallel = "no", ncpu = NULL, groups = NULL, distfn = dists_and_trans) {
  force(distfn)
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  ## if we get a constant outcome, then p=1.
  ## no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  thetreat <- fmla_vars[[2]]
  thedat <- copy(dat) #data.table::copy(dat) ## to avoid making changes outside the function
  ## stopifnot(is.data.table(thedat))
  outcome_names <- c(theresponse, "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "mhY", "rankY")
  # if(is.factor(thedat[[thetreat]])){
  #  thedat[,ZN:={tmp<-Rfast::as_integer(levels(get(thetreat)))[get(thetreat)];tmp-1}]
  # } else {
  #  thedat[,ZN:=get(thetreat)]
  #    }
  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := c(distfn(get(theresponse)), list(Rank(get(theresponse)))), by = get(theblock)]
    ## If one of the test statistics is constant, drop it.
    ## https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    # anyconstant_cols <-  whichAreConstant(thedat[,.SD,.SDcols=outcome_names], verbose=FALSE)
    # if(length(anyconstant_cols)>0) { outcome_names <- outcome_names[-anyconstant_cols] }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    ## This next is faster than doing it in two lines
    thedat[, outcome_names[-1] := c(distfn(get(theresponse)), list(Rank(get(theresponse))))]
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }
  newfmla <- as.formula(newfmla_text)
  if (is.null(simthresh) | nrow(thedat) > simthresh) {
    thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "quadratic"))[[1]]
    ## btfmla <- paste(thetreat,"~",paste(outcome_names,collapse="+"),"+strata(",theblock,")",sep="")
    ## thep2 <- RItools::balanceTest(as.formula(btfmla),data=thedat,report="chisquare.test")$overall[theblock,"p.value"]
    ## microbenchmark(pvalue(independence_test(newfmla, data = thedat, teststat = "quadratic"))[[1]],
    ##               RItools::balanceTest(as.formula(btfmla),data=thedat,report="chisquare.test")$overall[theblock,"p.value"],
    ##               times=20)
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- parallel::detectCores()
    }
    if (parallel == "no") {
      ncpu <- 1
    }
    thep <- pvalue(independence_test(newfmla,
      data = thedat,
      teststat = "quadratic",
      distribution = approximate(
        nresample = sims,
        parallel = parallel,
        ncpus = ncpu
      )
    ))[[1]]
  }
  return(as.numeric(thep))
}
