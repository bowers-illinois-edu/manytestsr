# Functions for splitting

#' Splitting function: K-Means Clustering
#'
#' A splitting function takes block ids and a block ordering
#' vector (or vectors) and produces a factor that assigns some block ids to one
#' group or another group.
#' @param bid Block id
#' @param x A vector that we can use to order the blocks
#' @return A factor categorizing blocks into groups.
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom Ckmeans.1d.dp Ckmeans.1d.dp
#' @export
splitCluster <- function(bid, x) {
  stopifnot("x must be numeric or integer" = is.numeric(x))

  if (length(unique(x)) == 1) {
    # Random splits used with stop_splitby_constant=FALSE otherwise find_blocks should stop before this
    group <- factor(sample(rep_len(c(0, 1), length.out = length(x))))
    return(group)
  }
  if (length(x) == 2) {
    # If there are only two values for x, split into two groups.
    group <- factor(c(0, 1))
    return(group)
  }
  ## if (length(x) == 3) {
  ##  # if we only have 3 values for x, make two groups with 1 group with the highest value
  ##  # and the other group with the 2 other values.
  ##  mnx <- fastMean(x)
  ##  rank_dists <- rank(abs(x - mnx))
  ##  group <- factor(as.numeric(rank_dists == max(rank_dists)))
  ##  return(group)
  ## }

  ### Trying to handle some edge
  ### cases with this function. Mostly KMeans_cpp works well but kmeans handles
  ### some cases where KMeans_cpp throws an error
  ## clus <- tryCatch(KMeans_rcpp(as.matrix(x), clusters = 2, num_init = 2,
  ##        initializer = "optimal_init")$clusters, error = function(e) {
  ##    kmeans(x, centers = 2)$cluster})

  ## Approach 2:
  clus <- Ckmeans.1d.dp(x, k = 2)$cluster
  group <- factor(as.numeric(clus == 1))
  # names(group) <- bid (no longer necessary)
  return(group)
}

#' Splitting function: Approx Equal Splits
#'
#' A splitting function takes block ids and block ordering
#' vector (or vectors) and produces a factor that assigns some block ids to one
#' group or another group such that the sum of x within one group is approximately equal to the sum of x in the other group. If x is something like block-size then this equalized the total observations in each split. If x is a covariate, then this splitting approach may make less conceptual sense --- if x <- c(1,2,3,4) it would put 1 and 3 in one group and 2 and 4 in another group.
#' @param bid Block id
#' @param x A vector that we can use to order the blocks
#' @return A factor categorizing blocks into groups.
#' @export
splitEqualApprox <- function(bid, x) {
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  names(x) <- bid
  sortedx <- sort(x)
  group1 <- sort(x)[seq_len(length(x)) %% 2 == 0]
  group <- factor(as.numeric(bid %in% names(group1)))
  return(group)
}


#' Splitting function: Leave One Out
#'
#' A splitting function takes block ids and block ordering
#' vector (or vectors) and produces a factor that assigns some block ids to one
#' group or another group.
#' @param bid Block id
#' @param x A vector that we can use to order the blocks. A block-level value.
#' Like N or harmonic mean weight of the block.
#' @return A factor categorizing blocks into groups.
#' @export
splitLOO <- function(bid, x) {
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  # require(data.table)
  # We only want to return one block and a time. So need to avoid ties.
  # Using random choice among equal ranks --- so this amounts to random splits when find_blocks doesn't stop the splitting because of stop_splitby_constant=TRUE
  frankx <- frank(x, ties.method = "random") # frank from data.table
  group <- factor(as.numeric(frankx < max(frankx)))
  # names(x) <- bid
  return(group)
}


#' A set of pre-specified splits
#'
#' This function does binary splits using a factor variable with dots separating the names of the subgroups. If there are more than two subgroups at any level, it makes one group from the largest subgroup and another group from the rest. If there are multiple subgroups the same size, it chooses the first subgroup by order of the levels of the factor.
#'
#' @param bid Block id
#' @param x Is a a factor with levels like "state.district.school". The splits will occur from left to right depending on whether there is existing variation at that level
#' @importFrom stringi stri_split_regex
#' @export
splitSpecifiedFactor <- function(bid, x) {
  stopifnot("x must be a factor" = is.factor(x))
  stopifnot("The factor must have categories separately by dots '.' and have at least one such separation." = stri_count_fixed(x, ".") > 0)
  if (length(unique(x)) == 1) {
    # Random splits used with stop_splitby_constant=FALSE otherwise find_blocks should stop before this
    group <- factor(sample(rep_len(c(0, 1), length.out = length(x))))
    return(group)
  }
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  x_split <- stri_split_regex(x, "\\.", simplify = TRUE)
  # Which  column varies (recalling that after previous splits some columns may have no variance).
  # try dataPreparation::whichAreConstant() or grab that function's C code in the future since it is much faster than below
  which_varies <- apply(x_split, 2, function(x) {
    length(unique(x))
  })
  # Split on the first column with variance from the beginning
  split_on <- which(which_varies > 1)[1]
  # Choose the largest subgroup to split from the rest or
  # make a random choice if they are the same size
  subgroup_sizes <- table(x_split[, split_on])
  ## Looks like which.max chooses either the max or the first of the max (in the order of the levels of the factor)
  which_group <- names(which.max(subgroup_sizes))
  # if (is.na(split_on)) {
  #  group <- factor(rep_len(0, length.out=length(x)))
  # } else {
  # as.numeric of factors creates a vector that starts a 1
  # group <- factor(as.numeric(x_split[, split_on] == x_split[1, split_on]))
  group <- factor(as.numeric(x_split[, split_on] == which_group))
  # }
  return(group)
}


#' A set of pre-specified splits
#'
#' This function allows for more than two splits at each level
#'
#' @param bid Block id
#' @param x Is a a factor with levels like "state.district.school". The splits will occur from left to right depending on whether there is existing variation at that level
#' @importFrom stringi stri_split_regex
#' @export
splitSpecifiedFactorMulti <- function(bid, x) {
  stopifnot(is.factor(x))
  stopifnot("The factor must have categories separately by dots '.' and have at least one such separation." = stri_count_fixed(x, ".") > 0)
  if (length(unique(x)) == 1) {
    # Random splits used with stop_splitby_constant=FALSE otherwise find_blocks should stop before this
    group <- factor(sample(rep_len(c(0, 1), length.out = length(x))))
    return(group)
  }
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  x_split <- stri_split_regex(x, "\\.", simplify = TRUE)
  # Which  column varies (recalling that after previous splits some columns may have no variance).
  # try dataPreparation::whichAreConstant() or grab that function's C code in the future since it is much faster than below
  which_varies <- apply(x_split, 2, function(x) {
    length(unique(x))
  })
  # Split on the first column with variance from the beginning
  split_on <- which(which_varies > 1)[1]
  # if (is.na(split_on)) {
  #  group <- factor(rep_len(0, length.out = length(x)))
  # } else {
  # as.numeric of factors creates a vector that starts at 1
  group <- factor(as.numeric(factor(x_split[, split_on])) - 1)
  # }
  return(group)
}


#' A set of pre-specified splits using a data.table object (Deprecate)
#'
#' @param bid Block id
#' @param x is a data.table object where each column from 1 to k further divides the blocks.
#' Column one should be highest level and each other column should be nested
#' @param bid is not used
#' @export
splitSpecified <- function(bid, x) {
  stopifnot("This splitting function requires a data.table object with more than one column" = is.data.table(x))
  stopifnot("Try a different splitting function if you have only a single attribute or criteria" = ncol(x) > 1)
  if (nrow(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  # Which  column varies (recalling that after previous splits some columns may have no variance).
  which_varies <- x[, lapply(.SD, function(x) {
    length(unique(x))
  }), .SDcols = names(x), drop = TRUE]
  # Split on the first column with variance from the beginning
  split_on <- names(which_varies)[as.vector(which_varies > 1)][1]
  # if (is.na(split_on)) {
  #  group <- factor(rep_len(0, length.out = nrow(x)))
  # } else {
  # as.numeric of factors creates a vector that starts a 1
  group <- as.numeric(factor(x[, get(split_on)])) - 1
  # }
  return(group)
}

## #### Deprecating the following function since it is slow and doesn't really add that much over splitEqualApprox
## #' Splitting function: Equal Splits
## #'
## #' A splitting function takes block ids and block ordering
## #' vector (or vectors) and produces a factor that assigns some block ids to one
## #' group or another group.
## #' @param bid Block id
## #' @param x A vector that we can use to order the blocks
## #' @return A factor categorizing blocks into groups.
## #' @importFrom nbpMatching nonbimatch distancematrix get.sets
## #' @export
## splitEqual <- function(bid, x) {
##   if (length(x) == 2) {
##     group <- factor(c(0, 1))
##     return(group)
##   }
##   names(x) <- bid
##   x <- sort(x)
##   oddx <- (length(x) %% 2) != 0
##   if (oddx) {
##     # if an odd number of blocks, exclude the least powerful one and add it later
##     smallx <- x[1] # x[x==min(x)][1]
##     x <- x[2:length(x)]
##   }
##   blockdists <- outer(x, x, function(x, y) {
##     abs(x - y)
##   })
##   blockdistsm <- distancematrix(blockdists)
##   sol <- nonbimatch(blockdistsm)
##   solsets <- get.sets(sol)
##   # Split into groups based on the pairs
##   group1sets <- names(sort(solsets)[(1:length(solsets) %% 2) == 0])
##   group <- factor(as.numeric((names(x) %in% group1sets)))
##   if (oddx) {
##     groupsum <- tapply(x, group, sum)
##     # the names of the groupnum object are the levels of the group factor
##     group <- factor(c(as.character(group), names(groupsum)[groupsum == min(groupsum)]),labels=c("0","1"))
##     # names(group) <- c(names(x),names(smallx))
##   }
##   return(group)
## }
##
# testSplitEqual <- splitEqual(trueBlockEffects$block,trueBlockEffects$hwt)
# stopifnot(length(testSplitEqual)==nrow(trueBlockEffects))
# stopifnot(sum(as.logical(testSplitEqual))==nrow(trueBlockEffects)/2)
# testSplitEqual <-
# splitEqual(trueBlockEffects$block[-1],trueBlockEffects$hwt[-1])
