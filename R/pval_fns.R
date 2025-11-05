## Functions for producing p-values / Test-statistics

#' P-value function: T-test
#'
#' These functions accept a data frame and perhaps test specific arguments
#' (like whether or not the test will be asymptotic or simulation based). It
#' produces a p-value.
#'
#' @param dat An object inheriting from class data.frame
#' @param fmla outcome~treatment factor | block factor (following `coin` API).

#' @param simthresh Below which number of total observations should the p-value
#' functions use permutations rather than asymptotic approximations

#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)

#' @param parallel Should the function use multicore processing for permutation
#' based testing. Default is no. But could be "snow" or "multicore" following
#' `approximate` in the coin package.

#' @param ncpu is the number of workers (for "snow") or cores (for "multicore").
#' @return A p-value
#' @examples
#' # Example using built-in data
#' data(example_dat, package = "manytestsr")
#'
#' # Test for treatment effect on Y1 within a single block
#' single_block <- subset(example_dat, blockF == "B080")
#' p_val <- pOneway(single_block, Y1 ~ trtF | blockF, parallel = "no")
#' print(p_val)
#'
#' # Test with permutation-based inference for small samples
#' p_val_perm <- pOneway(single_block, Y1 ~ trtF | blockF,
#'   simthresh = 100, sims = 500, parallel = "no"
#' )
#' print(p_val_perm)
#'
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

  if (is.null(simthresh) || nrow(dat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin:::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  thep <- pvalue(oneway_test(fmla, data = dat, distribution = thedist))[[1]]
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

#' @param simthresh Below which number of total observations should the p-value
#' functions use permutations rather than asymptotic approximations

#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)

#' @param parallel Should the function use multicore processing for permutation
#' based testing. Default is no. But could be "snow" or "multicore" following
#' `approximate` in the coin package.

#' @param ncpu is the number of workers (for "snow") or cores (for "multicore").
#' @return A p-value
#' @examples
#' # Example using Wilcoxon rank-sum test
#' data(example_dat, package = "manytestsr")
#'
#' # Test for treatment effect on Y1 within a single block
#' single_block <- subset(example_dat, blockF == "B080")
#' p_val <- pWilcox(single_block, Y1 ~ trtF | blockF, parallel = "no")
#' print(p_val)
#'
#' # Compare with permutation-based version
#' p_val_perm <- pWilcox(single_block, Y1 ~ trtF | blockF,
#'   simthresh = 100, sims = 500, parallel = "no"
#' )
#' print(p_val_perm)
#'
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

  if (is.null(simthresh) || nrow(dat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin:::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  thep <- pvalue(wilcox_test(fmla, data = dat, distribution = thedist))[[1]]

  return(as.numeric(thep))
}

#' P-value function: Independence Treatment Distance Test
#'
#' @description These functions accept a data frame and perhaps test specific
#' arguments (like whether or not the test will be asymptotic or simulation
#' based). It produces a p-value.
#'
#' @details For now, this  function  does an omnibus-style chi-square  test
#' using (1) the ratio  of  distances to controls to distances to treated
#' observations within block;  (2)  the  rank of  distances to  controls for
#' each unit; and (3) the raw outcome.
#'
#' Although the distances are calculated by block, our profiling suggests that
#' it is better to parallelize the distance creation `distfn` (done here in C++
#' in the `fastfns.cpp` file) rather than use the `data.table` approach of
#' `setDTthreads()`. So, here we assume that the threads for data.table are 1.

#' @param dat An object inheriting from class data.frame
#' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param simthresh is the size of the data below which we use direct permutations for p-values

#' @param distfn is  a function that produces one or more vectors (a data frame
#' or matrix) of the  same number of  rows as the dat

#' @param parallel is "no" then parallelization is not required, otherwise it
#' is "multicore" or "snow" in the call to `coin::independence_test()` (see
#' help for coin::approximate()). Also, if parallel is not "no" and
#' `adaptive_dist_function` is TRUE, then an openmp version of the distance
#' creation function is called using `ncpu` threads (or
#' `parallel::detectCores(logical=FALSE)` cores).

#' @param ncpu is number of cpus  to be used for parallel operation.

#' @return A p-value
#' @importFrom coin independence_test pvalue approximate exact asymptotic
#' @importFrom Rfast Rank
#' @importFrom parallel detectCores
#' @importFrom dataPreparation which_are_constant
#' @examples
#' \donttest{
#' # Example using distance-based independence test
#' data(example_dat, package = "manytestsr")
#' library(data.table)
#'
#' # Test for treatment effect using distance-based approach
#' single_block <- as.data.table(subset(example_dat, blockF == "B080"))
#' p_val <- pIndepDist(single_block, Y1 ~ trtF | blockF, parallel = "no")
#' print(p_val)
#'
#' # Test with different outcome variable
#' p_val2 <- pIndepDist(single_block, Y2 ~ trtF | blockF, parallel = "no")
#' print(p_val2)
#' }
#' @export
pIndepDist <- function(dat, fmla = YcontNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                       parallel = "yes", ncpu = NULL, distfn = fast_dists_and_trans_hybrid) {
  force(distfn)
  stopifnot(inherits(dat, "data.table"))
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  thetreat <- fmla_vars[[2]]
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

  thedat <- copy(dat)

  ### These must match the names of the functions used in src/dists_and_trans.cpp
  outcome_names <- c(theresponse, "mean_dist", "mean_rank_dist", "max_dist", "rankY", "tanhY")

  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := distfn(get(theresponse)), by = get(theblock)]
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
    thedat[, outcome_names[-1] := distfn(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }
  newfmla <- as.formula(newfmla_text)

  if (is.null(simthresh) || nrow(dat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin:::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  ## Quadratic combination is best when most/all of the test statistics move in
  ## the same direction. Gains power from combining. maxT is best when we think
  ## that perhaps only one of the set of scores will show an effect

  thep <- pvalue(independence_test(newfmla, data = thedat, teststat = "quadratic", distribution = thedist))[[1]]
  return(as.numeric(thep))
}

#' P-value function: Testing twice
#'
#' @description These functions accept a data frame and perhaps test specific
#' arguments (like whether or not the test will be asymptotic or simulation
#' based). It produces a p-value.
#'
#' @details For now, this  function  does an omnibus-style max-T test using (1)
#' the raw outcome and (2) a rank transformed raw outcome. Inspired by Rosenbaum (2008) on Testing Twicee
#'

#' @param dat An object inheriting from class data.frame
#' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param simthresh is the size of the data below which we use direct permutations for p-values

#' @param distfn is  a function that produces one or more vectors (a data frame
#' or matrix) of the  same number of  rows as the dat

#' @param parallel is "no" then parallelization is not required, otherwise it
#' is "multicore" or "snow" in the call to `coin::independence_test()` (see
#' help for coin::approximate()). Also, if parallel is not "no" and
#' `adaptive_dist_function` is TRUE, then an openmp version of the distance
#' creation function is called using `ncpu` threads (or
#' `parallel::detectCores(logical=FALSE)` cores).

#' @param ncpu is number of cpus  to be used for parallel operation.

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
  thetreat <- fmla_vars[[2]]
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

  thedat <- copy(dat)
  outcome_names <- c(theresponse, "rankY")
  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, rankY := frank(get(theresponse)), by = get(theblock)]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    thedat[, rankY := frank(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }
  newfmla <- as.formula(newfmla_text)
  if (is.null(simthresh) || nrow(dat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin:::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  thep <- pvalue(independence_test(newfmla, data = thedat, distribution = thedist, teststat = "quadratic"))[[1]]
  return(as.numeric(thep))
}

#' P-value function: Cauchy Combined Indepence Test
#'

#' @details This  function combines the p-values from univariate tests using
#' using (1)  the raw outcome, (2) the rank-transformed outcome, (3) the tanh
#' transformed raw outcome (another test statistic that is more sensitive when
#' there is skew), (4) the mean difference in  raw outcome Euclidean distances
#' bwtween treated and control observations within block;  (5) the mean
#' difference in ranked outcome Euclidean distances bwtween treated and control
#' observations; (6) a Hotelling-T style quadratic combination of the preceding
#' test statistics. Inspired by Rizzo and Székely's work on Euclidean distance
#' based testing and by Liu and Xie (2020) on the Cauchy Combination Test and
#' Hansen and Bowers (2008) on omnibus tests. Distance and ranks and other
#' transformations are all calculated by block when block is supplied in the
#' formula.

#' @param dat An object inheriting from class data.frame
#' @param fmla  A formula  appropriate to the function. Here it should  be something like outcome~treatment|block
#' @param sims Either NULL (meaning use an asymptotic reference dist) or a
#' number (meaning sampling from the randomization distribution implied by the
#' formula)
#' @param simthresh is the size of the data below which we use direct permutations for p-values

#' @param distfn is  a function that produces one or more vectors (a data frame
#' or matrix) of the  same number of  rows as the dat.  Here until further
#' notice users should leave it at `mean_dist_raw_rank_tanh` since that function is
#' purpose built for this function.

#' @param parallel is "no" then parallelization is not required, otherwise it
#' is "multicore" or "snow" in the call to `coin::independence_test()` (see
#' help for coin::approximate()). Also, if parallel is not "no" and
#' `adaptive_dist_function` is TRUE, then an openmp version of the distance
#' creation function is called using `ncpu` threads (or
#' `parallel::detectCores(logical=FALSE)` cores).

#' @param ncpu is number of cpus  to be used for parallel operation.
#' @return A p-value
#' @importFrom coin independence_test pvalue approximate exact asymptotic oneway_test wilcox_test
#' @importFrom Rfast Rank
#' @importFrom parallel detectCores
#' @importFrom dataPreparation which_are_constant
#' @export
pCombCauchyDist <- function(dat, fmla = YcontNorm ~ trtF | blockF, simthresh = 20, sims = 1000,
                            parallel = "no", ncpu = NULL, distfn = fast_dists_and_trans_nomax_hybrid) {
  force(distfn)
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  thetreat <- fmla_vars[[2]]
  # if we get a constant outcome, then p=1.
  # no evidence against the null of no effects.
  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }

  if (is.null(ncpu) && parallel != "no") {
    ## This is used for both distfn but also for independence_test for it to do
    ## permutation testing in parallel in small blocks/samples

    ncpu <- detectCores()
  }

  thedat <- copy(dat)

  outcome_names <- c(theresponse, "mean_dist", "mean_rank_dist", "rankY", "tanhY")

  ## To enable the use of this function without a `|block` --- basically to
  ## make it easier to compare with the bottom-up approaches where we just test
  ## within each block. We do the following: (TODO)

  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := distfn(get(theresponse)), by = get(theblock)]
    # If one of the test statistics is constant, drop it.
    # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    ## Make a list of formulas that we will use to get the p-values for combination below
    the_fmlas <- lapply(outcome_names, function(ynm) {
      return(paste(ynm, "~", thetreat, "|", theblock, sep = ""))
    })
    the_fmlas[[5]] <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    ### If there is no block in the formula and this is a test that will be
    ### done within each block and not aggregated across them
    thedat[, outcome_names[-1] := distfn(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    the_fmlas <- lapply(outcome_names, function(ynm) {
      return(paste(ynm, "~", thetreat, sep = ""))
    })
    the_fmlas[[5]] <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  }

  ## Add a combined version in case the signal is only detectable this way

  if (is.null(simthresh) || nrow(thedat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin:::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  p_vals <- sapply(the_fmlas, function(thefmla) {
    p_tmp <- pvalue(independence_test(as.formula(thefmla), data = thedat, teststat = "quadratic", distribution = thedist))[[1]]
    return(p_tmp)
  })

  ## Now we will calculate separate p-values and combine them using the Cauchy combination method (Liu and Xie 2019)
  ## Cauchy Combination Test: A Powerful Test With Analytic p-Value Calculation Under Arbitrary Dependency Structures
  ## This function goes kind of crazy at 0 and 1 so add or substract a tiny amount from either side if needed
  acat_pvalue <- function(p_values) {
    if (any(p_values == 0 | p_values == 1, na.rm = TRUE)) {
      #      stop("Input p-values must lie strictly inside (0,1).")
      p_values[p_values == 0] <- p_values[p_values == 0] + .Machine$double.eps
      p_values[p_values == 1] <- p_values[p_values == 1] - .Machine$double.eps
    }
    ## pi is pi FYI
    ## Since this procedure creates weird results when p is very near 1 or p very near 0, we create
    ## the two sided p-values using one-sided tests and then double the minimum.
    p_1 <- p_values / 2
    t_stat <- mean(tan((0.5 - p_1) * pi)) # equal weights
    upper_p <- 0.5 - atan(t_stat) / pi # final ACAT p-value
    ## We want a 2 sided test
    return(2 * min(c(upper_p, 1 - upper_p)))
  }

  thep <- acat_pvalue(p_vals)

  return(as.numeric(thep))
}
