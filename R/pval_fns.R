# Functions for producing p-values / Test-statistics

#' P-value function: T-test
#'
#' These functions accept a data frame and perhaps test specific arguments
#' (like whether or not the test will be asymptotic or simulation based). It
#' produces a p-value.
#'
#' @param dat An object inheriting from class data.frame
#' @param fmla outcome~treatment factor | block factor (following `coin` API).
#' @param simthresh Below which number of total observations should the p-value functions use permutations rather than asymptotic approximations
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param parallel Should the function use multicore processing for permutation based testing. Default is no. But could be "snow" or "multicore" following `approximate` in the coin package.
#' @param ncpu is the number of workers (for "snow") or cores (for "multicore").
#' @return A p-value
#' @importFrom coin oneway_test pvalue approximate exact asymptotic
#' @importFrom parallel detectCores
#' @export
pOneway <- function(dat, fmla = YContNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                    parallel = "no", ncpu = NULL) {
  theresponse <- all.vars(fmla)[attr(terms(fmla), "response")]
  # if we get a constant outcome, then p=1.
  # no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  if (is.null(simthresh) | nrow(dat) > simthresh) {
    thep <- pvalue(oneway_test(fmla, data = dat))[[1]]
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- detectCores()
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

#' P-value function: Wilcox Test
#'
#' These functions accept a data frame and perhaps test specific arguments
#' (like whether or not the test will be asympotic or simulation based). It
#' produces a p-value.
#'
#' @param dat An object inheriting from class data.frame
#' @param fmla outcome~treatment factor | block factor (following `coin` API).
#' @param simthresh Below which number of total observations should the p-value functions use permutations rather than asymptotic approximations
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param parallel Should the function use multicore processing for permutation based testing. Default is no. But could be "snow" or "multicore" following `approximate` in the coin package.
#' @param ncpu is the number of workers (for "snow") or cores (for "multicore").
#' @return A p-value
#' @importFrom coin wilcox_test pvalue approximate exact asymptotic
#' @importFrom parallel detectCores
#' @export
pWilcox <- function(dat, fmla = YContNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                    parallel = "no", ncpu = NULL) {
  theresponse <- all.vars(fmla)[attr(terms(fmla), "response")]
  # if we get a constant outcome, then p=1.
  # no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }
  if (is.null(simthresh) | nrow(dat) > simthresh) {
    thep <- pvalue(wilcox_test(fmla, data = dat))[[1]]
  } else {
    if (is.null(ncpu) & parallel != "no") {
      ncpu <- detectCores()
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

#' P-value function: Independence Treatment Distance Test
#'
#' @description
#' These functions accept a data frame and perhaps test specific arguments
#' (like whether or not the test will be asymptotic or simulation based). It
#' produces a p-value.
#'
#' @details
#' For now, this  function  does an omnibus-style chi-square  test using (1)
#' the ratio  of  distances to controls to distances to treated observations
#' within block;  (2)  the  rank of  distances to  controls for each unit; and
#' (3) the raw outcome.
#'
#' Although the distances are calculated by block, our profiling suggests that it is better to parallelize the distance creation `distfn` (done here in C++ in the `fastfns.cpp` file) rather than use the `data.table` approach of `setDTthreads()`. So, here we assume that the threads for data.table are 1.
#' @param dat An object inheriting from class data.frame
#' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param simthresh is the size of the data below which we use direct permutations for p-values
#' @param groups is a vector defining the groups within which the inter-unit  distances are calculated. Not used here.
#' @param distfn is  a function that produces one or more vectors (a data frame or matrix) of the  same number of  rows as the dat
#' @param parallel is "no" then parallelization is not required, otherwise it is "multicore" or "snow" in the call to `coin::independence_test()` (see  help for coin::approximate()). Also, if parallel is not "no" and `adaptive_dist_function` is TRUE, then an openmp version of the distance creation function is called using `ncpu` threads (or `parallel::detectCores(logical=FALSE)` cores).
#' @param ncpu is number of cpus  to be used for parallel operation.
#' @param adaptive_dist_function is TRUE if the distance calculation function should be chosen using previous benchmarks. See the code.
#' @return A p-value
#' @importFrom coin independence_test pvalue approximate exact asymptotic
#' @importFrom Rfast Rank
#' @importFrom parallel detectCores
#' @importFrom dataPreparation which_are_constant
#' @export
pIndepDist <- function(dat, fmla = YcontNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                       parallel = "yes", ncpu = NULL, groups = NULL, distfn = dists_and_trans, adaptive_dist_function = TRUE) {
  force(distfn)
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  # if we get a constant outcome, then p=1.
  # no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }

  dat_size <- nrow(dat)

  if (is.null(ncpu) & parallel != "no") {
    ## This is used for both distfn but also for independence_test for it to do permutation testing in parallel in small blocks/samples
    ncpu <- detectCores()
  }

  ## the dists_and_trans function is very fast but very memory intensive. When n is large it can require GB versus KB
  ## This should be about the sizes of the blocks not the overall size. TODO
  if (adaptive_dist_function) {
    if (dat_size <= 100) {
      distfn <- dists_and_trans
    }
    if (parallel != "no" & dat_size > 100) {
      distfn <- fast_dists_and_trans_by_unit_arma_parR(threads = ncpu)
    }
    if (parallel == "no" & dat_size > 100) {
      distfn <- fast_dists_and_trans_by_unit_arma2
    }
  }

  thetreat <- fmla_vars[[2]]
  thedat <- copy(dat)
  #  outcome_names <- c(theresponse, "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "mhdist", "rankx", "mnsqrtdist", "hubmn", "tanhx")
  # outcome_names <- c(theresponse, "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")
  outcome_names <- c(theresponse, "mndist", "mndistRank0", "maxdist", "rankY", "tanhx")
  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := distfn(get(theresponse), Z = 1), by = get(theblock)]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    # This next is faster than doing it in two lines
    thedat[, outcome_names[-1] := distfn(get(theresponse), Z = 1)]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }
  newfmla <- as.formula(newfmla_text)
  if (is.null(simthresh) | nrow(thedat) > simthresh) {
    suppressWarnings(
      # thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "maximum", alternative="greater"))[[1]]
      thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "quadratic"))[[1]]
    )
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    suppressWarnings(
      thep <- pvalue(independence_test(newfmla,
        data = thedat,
        teststat = "quadratic",
        ## 	alternative= "greater",
        distribution = approximate(
          nresample = sims,
          parallel = parallel,
          ncpus = ncpu
        )
      ))[[1]]
    )
  }
  return(as.numeric(thep))
}

#' P-value function: Testing twice
#'
#' @description
#' These functions accept a data frame and perhaps test specific arguments
#' (like whether or not the test will be asymptotic or simulation based). It
#' produces a p-value.
#'
#' @details
#' For now, this  function  does an omnibus-style max-T test using (1) the raw outcome and (2) a rank transformed raw outcome.
#'
#' @param dat An object inheriting from class data.frame
#' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param simthresh is the size of the data below which we use direct permutations for p-values
#' @param groups is a vector defining the groups within which the inter-unit  distances are calculated. Not used here.
#' @param distfn is  a function that produces one or more vectors (a data frame or matrix) of the  same number of  rows as the dat
#' @param parallel is "no" then parallelization is not required, otherwise it is "multicore" or "snow" in the call to `coin::independence_test()` (see  help for coin::approximate()). Also, if parallel is not "no" and `adaptive_dist_function` is TRUE, then an openmp version of the distance creation function is called using `ncpu` threads (or `parallel::detectCores(logical=FALSE)` cores).
#' @param ncpu is number of cpus  to be used for parallel operation.
#' @param adaptive_dist_function is TRUE if the distance calculation function should be chosen using previous benchmarks. See the code.
#' @return A p-value
#' @importFrom coin independence_test pvalue approximate exact asymptotic
#' @importFrom Rfast Rank
#' @importFrom parallel detectCores
#' @importFrom dataPreparation which_are_constant
#' @export
pTestTwice <- function(dat, fmla = YcontNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                       parallel = "yes", ncpu = NULL, groups = NULL) {
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  # if we get a constant outcome, then p=1.
  # no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }

  dat_size <- nrow(dat)

  if (is.null(ncpu) && parallel != "no") {
    ## This is used for both distfn but also for independence_test for it to do
    ## permutation testing in parallel in small blocks/samples

    ncpu <- detectCores()
  }

  thetreat <- fmla_vars[[2]]
  thedat <- copy(dat)
  outcome_names <- c(theresponse, "rankY")
  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := frank(get(theresponse)), by = get(theblock)]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    # This next is faster than doing it in two lines
    thedat[, outcome_names[-1] := frank(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }
  newfmla <- as.formula(newfmla_text)
  if (is.null(simthresh) | nrow(thedat) > simthresh) {
    suppressWarnings(
      # thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "maximum", alternative="greater"))[[1]]
      thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "maximum"))[[1]]
    )
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    suppressWarnings(
      thep <- pvalue(independence_test(newfmla,
        data = thedat,
        teststat = "maximum",
        ## 	alternative= "greater",
        distribution = approximate(
          nresample = sims,
          parallel = parallel,
          ncpus = ncpu
        )
      ))[[1]]
    )
  }
  return(as.numeric(thep))
}
