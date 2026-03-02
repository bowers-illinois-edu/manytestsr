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
    thedist <- coin::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
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
    thedist <- coin::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
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
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    # This next is faster than doing it in two lines
    thedat[, outcome_names[-1] := distfn(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
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
    thedist <- coin::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
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
#' @param parallel is "no" then parallelization is not required, otherwise it
#' is "multicore" or "snow" in the call to `coin::independence_test()` (see
#' help for coin::approximate()). Also, if parallel is not "no" and
#' `adaptive_dist_function` is TRUE, then an openmp version of the distance
#' creation function is called using `ncpu` threads (or
#' `parallel::detectCores(logical=FALSE)` cores).
#' @param ncpu is number of cpus  to be used for parallel operation.
#' @param groups Currently unused parameter, reserved for future functionality
#'
#' @return A p-value
#' @importFrom coin independence_test pvalue approximate exact asymptotic
#' @importFrom Rfast Rank
#' @importFrom parallel detectCores
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
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
    if (length(anyconstant_cols) > 0) {
      outcome_names <- outcome_names[-anyconstant_cols]
    }
    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  } else {
    theblock <- NULL
    thedat[, rankY := frank(get(theresponse))]
    # If one of the test statistics is constant, drop it.
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
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
    thedist <- coin::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
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
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
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
    anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
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
    thedist <- coin::approximate(object, nresample = sims, parallel = parallel, ncpus = ncpu)
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


#' P-value function: Polynomial Rank Score Test
#'
#' @description Tests Fisher's sharp null of no treatment effects using
#' multiple polynomial rank score functions simultaneously via
#' \code{coin::independence_test()}. This provides adaptive sensitivity
#' to effects at different parts of the outcome distribution without
#' pre-committing to a single scoring.
#'
#' @details For each value of r in \code{r_vec}, the function computes
#' within-block polynomial rank scores: \code{(rank(Y) / (n_b + 1))^(r-1)}
#' where \code{n_b} is the block size. These scores are passed as a
#' multivariate response to \code{coin::independence_test()}, which
#' computes a joint test statistic that accounts for the correlation
#' structure among the score functions.
#'
#' With r = 2 the scores are nearly linear in rank (Wilcoxon-like).
#' Larger r values place increasing weight on observations with high
#' ranks, providing sensitivity to treatment effects concentrated in the
#' upper tail of the outcome distribution.
#'
#' This function tests the same sharp null as \code{\link{pOneway}} and
#' \code{\link{pIndepDist}} but uses polynomial rank scores from the
#' combined Stephenson rank test framework (Kim, Li, and Bowers).
#' For quantile-of-effects hypotheses (whether the k-th largest
#' individual effect exceeds a threshold), see
#' \code{\link{pCombStephenson}}.
#'
#' @param dat An object inheriting from class data.frame (will be
#'   coerced to data.table internally if needed).
#' @param fmla A formula: outcome ~ treatment | block.
#' @param r_vec Numeric vector. Polynomial rank score parameters.
#'   Default \code{c(2, 6, 10)}. r = 2 gives Wilcoxon-like scores;
#'   larger r emphasizes extreme ranks.
#' @param teststat Character. Type of test statistic passed to
#'   \code{coin::independence_test()}. \code{"quadratic"} (default)
#'   gives an omnibus chi-squared statistic that pools across score
#'   functions; \code{"maximum"} takes the max of the standardized
#'   statistics across score functions.
#' @param simthresh Integer. Below this number of observations, use
#'   permutation-based inference instead of asymptotic approximation.
#' @param sims Integer. Number of permutations when using simulation-based
#'   inference.
#' @param parallel Character. \code{"no"} (default), \code{"snow"}, or
#'   \code{"multicore"}, passed to \code{coin::approximate()}.
#' @param ncpu Integer. Number of CPUs for parallel operation.
#'
#' @return A p-value (numeric scalar).
#'
#' @examples
#' # Example using built-in data
#' data(example_dat, package = "manytestsr")
#' library(data.table)
#'
#' # Test for treatment effect using polynomial rank scores
#' single_block <- as.data.table(subset(example_dat, blockF == "B080"))
#' p_val <- pPolyRank(single_block, Y1 ~ trtF | blockF)
#' print(p_val)
#'
#' @importFrom coin independence_test pvalue approximate asymptotic
#' @export
pPolyRank <- function(dat, fmla = Y ~ trtF | blockF,
                      r_vec = c(2, 6, 10),
                      teststat = "quadratic",
                      simthresh = 20, sims = 1000,
                      parallel = "no", ncpu = NULL) {
  stopifnot(inherits(dat, "data.frame"))

  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  thetreat <- fmla_vars[[2]]

  if (length(unique(dat[[theresponse]])) < 2) {
    return(1)
  }

  thedat <- data.table::copy(data.table::as.data.table(dat))

  # Compute polynomial rank scores within blocks for each r
  score_names <- paste0("poly_r", r_vec)
  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    for (i in seq_along(r_vec)) {
      r <- r_vec[i]
      thedat[, (score_names[i]) := {
        rk <- rank(get(theresponse), ties.method = "average")
        nb <- .N
        (rk / (nb + 1))^(r - 1)
      }, by = get(theblock)]
    }
  } else {
    theblock <- NULL
    for (i in seq_along(r_vec)) {
      r <- r_vec[i]
      thedat[, (score_names[i]) := {
        rk <- rank(get(theresponse), ties.method = "average")
        nb <- .N
        (rk / (nb + 1))^(r - 1)
      }]
    }
  }

  # Drop constant score columns
  anyconstant_cols <- which_are_constant(thedat[, .SD, .SDcols = score_names], verbose = FALSE)
  if (length(anyconstant_cols) > 0) {
    score_names <- score_names[-anyconstant_cols]
  }
  if (length(score_names) == 0) {
    return(1)
  }

  # Build formula: poly_r2 + poly_r6 + poly_r10 ~ treatment | block
  if (!is.null(theblock)) {
    newfmla_text <- paste(paste(score_names, collapse = "+"), "~", thetreat, "|", theblock)
  } else {
    newfmla_text <- paste(paste(score_names, collapse = "+"), "~", thetreat)
  }
  newfmla <- as.formula(newfmla_text)

  if (is.null(simthresh) || nrow(dat) > simthresh) {
    thedist <- "asymptotic"
  } else {
    if (parallel == "no") {
      ncpu <- 1
    }
    thedist <- coin::approximate(nresample = sims, parallel = parallel, ncpus = ncpu)
  }

  thep <- pvalue(independence_test(newfmla, data = thedat,
    teststat = teststat, distribution = thedist))[[1]]
  return(as.numeric(thep))
}


#' P-value function: Combined Stephenson Rank Test
#'
#' @description Wraps \code{CMRSS::pval_comb_block()} to provide a
#' consistent formula interface for the combined Stephenson rank test
#' of Kim, Li, and Bowers. Tests quantile-of-effects hypotheses: whether
#' the k-th largest individual treatment effect exceeds a threshold c.
#'
#' @details The combined test uses polynomial rank scores at multiple
#' tuning parameters (controlled by \code{r_vec}) and takes the maximum
#' of the standardized statistics, avoiding the need to choose a single
#' tuning parameter. This wrapper builds the \code{methods.list.all}
#' configuration that \code{CMRSS::pval_comb_block()} expects.
#'
#' The \code{k} parameter indexes the k-th largest treatment effect,
#' where tau_(1) >= tau_(2) >= ... >= tau_(n). The test asks whether
#' tau_(k) exceeds \code{c}. Internally, CMRSS sets the top
#' \code{min(m, n-k)} treated units' adjusted outcomes to infinity.
#' The test is non-trivial only when \code{k > n - m} (where n is
#' total units and m is number treated), because otherwise all
#' treated units receive infinity and the test statistic is constant
#' across all permutations.
#'
#' The default \code{k = n} and \code{c = 0} tests Fisher's sharp
#' null of no effects for any unit. At k = n no treated units receive
#' infinity, so the test statistic reduces to a standard stratified
#' rank-sum statistic — but combined across multiple polynomial score
#' functions (controlled by \code{r_vec}), which provides adaptive
#' sensitivity without pre-committing to a single scoring. Smaller k
#' values test quantile-of-effects hypotheses: whether at least
#' \code{n - k + 1} units have effects exceeding \code{c}.
#'
#' @param dat An object inheriting from class data.frame (will be
#'   coerced to data.table internally if needed).
#' @param fmla A formula: outcome ~ treatment | block.
#' @param k Integer. The quantile index: test whether the k-th largest
#'   individual effect exceeds \code{c}. Default is \code{n} (the
#'   total number of units), which tests Fisher's sharp null of no
#'   effects. Must satisfy \code{k > n - m} for the test to be
#'   non-trivial. Smaller k tests quantile-of-effects hypotheses:
#'   whether at least \code{n - k + 1} units have effects exceeding
#'   \code{c}.
#' @param c Numeric. The threshold for the null hypothesis. Default 0.
#' @param r_vec Integer vector. Polynomial rank score tuning parameters.
#'   Default \code{c(2, 6, 10)}. These control the weight placed on
#'   different parts of the rank distribution: r = 2 is close to
#'   Wilcoxon, larger r emphasizes extreme ranks.
#' @param weight_name Character. Weighting scheme across blocks.
#'   \code{"asymp.opt"} (default) uses asymptotically optimal weights;
#'   \code{"dis.free"} uses distribution-free weights.
#' @param null_max Integer. Number of permutations for the
#'   randomization null distribution. Default 100000.
#' @param opt_method Character. Optimization method for the test
#'   statistic. Default \code{"ILP_auto"} (integer linear programming
#'   with automatic solver selection). See \code{?CMRSS::pval_comb_block}
#'   for options.
#' @param comb_method Integer. 1 = aggregate across strata then combine
#'   across methods (default); 2 = combine across methods within each
#'   stratum then aggregate (often more powerful).
#'
#' @return A p-value (numeric scalar).
#'
#' @examples
#' \donttest{
#' # Requires CMRSS package:
#' # remotes::install_github("bowers-illinois-edu/CMRSS")
#' if (requireNamespace("CMRSS", quietly = TRUE)) {
#'   data(example_dat, package = "manytestsr")
#'   library(data.table)
#'   idat <- as.data.table(example_dat)
#'   p <- pCombStephenson(idat, Y1 ~ trtF | blockF)
#'   print(p)
#' }
#' }
#'
#' @export
pCombStephenson <- function(dat, fmla = Y ~ trtF | blockF,
                            k = NULL, c = 0,
                            r_vec = c(2, 6, 10),
                            weight_name = "asymp.opt",
                            null_max = 10^5,
                            opt_method = "ILP_auto",
                            comb_method = 1) {
  if (!requireNamespace("CMRSS", quietly = TRUE)) {
    stop(
      "Package 'CMRSS' is required for this function. ",
      "Install it with: remotes::install_github('bowers-illinois-edu/CMRSS')",
      call. = FALSE
    )
  }

  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  thetreat <- fmla_vars[[2]]

  # Block is required — the combined Stephenson test aggregates across strata
  if (length(fmla_vars) < 3) {
    stop("pCombStephenson requires a block variable in the formula ",
         "(e.g., Y ~ trt | block).", call. = FALSE)
  }
  theblock <- fmla_vars[[3]]

  Y <- dat[[theresponse]]
  Z <- as.numeric(dat[[thetreat]])
  block <- factor(dat[[theblock]])

  # If treatment is a factor with levels like "0"/"1", convert properly
  if (is.factor(dat[[thetreat]])) {
    Z <- as.numeric(levels(dat[[thetreat]])[dat[[thetreat]]])
  }
  stopifnot(all(Z %in% c(0, 1)))

  n <- length(Z)
  m <- sum(Z)

  # Default k = n: test Fisher's sharp null of no effects.
  # At k = n, no treated units receive xi = Inf, and the test statistic
  # reduces to a standard stratified rank-sum (combined across score
  # functions in r_vec). For quantile-of-effects hypotheses, the user
  # can set k < n; k must exceed n - m for a non-trivial test.
  if (is.null(k)) {
    k <- n
  }
  if (k <= n - m) {
    warning(
      "k = ", k, " is at most n - m = ", n - m,
      ". The test statistic will be degenerate (all treated units ",
      "receive xi = Inf). Set k > ", n - m, " for a non-trivial test.",
      call. = FALSE
    )
  }

  # Build methods.list.all: H tuning parameters x B blocks
  # Each entry specifies a polynomial rank score configuration
  B <- nlevels(block)
  methods.list.all <- lapply(r_vec, function(r) {
    lapply(seq_len(B), function(b) {
      list(name = "Polynomial", r = r, std = TRUE, scale = FALSE)
    })
  })

  result <- CMRSS::pval_comb_block(
    Z = Z, Y = Y, k = k, c = c,
    block = block,
    methods.list.all = methods.list.all,
    weight.name = weight_name,
    null.max = null_max,
    statistic = FALSE,
    opt.method = opt_method,
    comb.method = comb_method
  )

  return(as.numeric(result))
}
