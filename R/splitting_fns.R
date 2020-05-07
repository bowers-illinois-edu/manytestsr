## Functions for splitting

##' Splitting function: K-Means Clustering
##'
##' A splitting function takes block ids and a block ordering
##' vector (or vectors) and produces a factor that assigns some block ids to one
##' group or another group.
##' @param bid Block id
##' @param x A vector that we can use to order the blocks
##' @return A factor categorizing blocks into groups.
##' @export
splitCluster <- function(bid, x) {
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  clus <- kmeans(x, centers = 2, nstart = 10)
  group <- factor(as.numeric(clus$cluster == 1))
  names(group) <- bid
  return(group)
}
## debug(splitEqualApprox) testSplitEqualApprox <-
## splitEqualApprox(bdat$block,bdat$hwt)



##' Splitting function: Equal Splits
##'
##' A splitting function takes block ids and block ordering
##' vector (or vectors) and produces a factor that assigns some block ids to one
##' group or another group.
##' @param bid Block id
##' @param x A vector that we can use to order the blocks
##' @return A factor categorizing blocks into groups.
##' @export
splitEqual <- function(bid, x) {
  require(nbpMatching, quietly = TRUE)
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  names(x) <- bid
  x <- sort(x)
  oddx <- (length(x) %% 2) != 0
  if (oddx) {
    ## if an odd number of blocks, exclude the least powerful one and add it later
    smallx <- x[1] ## x[x==min(x)][1]
    x <- x[2:length(x)]
  }
  blockdists <- outer(x, x, function(x, y) {
    abs(x - y)
  })
  blockdistsm <- nbpMatching::distancematrix(blockdists)
  sol <- nbpMatching::nonbimatch(blockdistsm)
  solsets <- nbpMatching::get.sets(sol)
  ## Split into groups based on the pairs
  group1sets <- names(sort(solsets)[(1:length(solsets) %% 2) == 0])
  group <- factor(as.numeric((names(x) %in% group1sets)))
  if (oddx) {
    groupsum <- tapply(x, group, sum)
    ## the names of the groupnum object are the levels of the group factor
    group <- factor(c(as.character(group), names(groupsum)[groupsum == min(groupsum)]))
    ## names(group) <- c(names(x),names(smallx))
  }
  return(group)
}

## testSplitEqual <- splitEqual(trueBlockEffects$block,trueBlockEffects$hwt)
## stopifnot(length(testSplitEqual)==nrow(trueBlockEffects))
## stopifnot(sum(as.logical(testSplitEqual))==nrow(trueBlockEffects)/2)
## testSplitEqual <-
## splitEqual(trueBlockEffects$block[-1],trueBlockEffects$hwt[-1])


##' Splitting function: Approx Equal Splits
##'
##' A splitting function takes block ids and block ordering
##' vector (or vectors) and produces a factor that assigns some block ids to one
##' group or another group.
##' @param bid Block id
##' @param x A vector that we can use to order the blocks
##' @return A factor categorizing blocks into groups.
##' @export
splitEqualApprox <- function(bid, x) {
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  names(x) <- bid
  sortedx <- sort(x)
  group1 <- sort(x)[1:length(x) %% 2 == 0]
  group <- factor(as.numeric(bid %in% names(group1)))
  return(group)
}
## debug(splitEqualApprox) testSplitEqualApprox <-
## splitEqualApprox(bdat$block,bdat$hwt)


##' Splitting function: Leave One Out
##'
##' A splitting function takes block ids and block ordering
##' vector (or vectors) and produces a factor that assigns some block ids to one
##' group or another group.
##' @param bid Block id
##' @param x A vector that we can use to order the blocks. A block-level value.
##' Like N or harmonic mean weight of the block.
##' @return A factor categorizing blocks into groups.
##' @export
splitLOO <- function(bid, x) {
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  ## require(data.table)
  ## We only want to return one block and a time. So need to avoid ties.
  frankx <- frank(x, ties.method = "first") ## frank from data.table
  group <- factor(as.numeric(frankx < max(frankx)))
  ## names(x) <- bid
  return(group)
}


##' A set of pre-specified splits
##'
##' @param x Is a a factor with levels like "state.district.school". The splits will occur from left to right depending on whether there is existing variation at that level
##' @importFrom stringi stri_split_regex
##' @export
splitSpecifiedFactor <- function(bid, x) {
  stopifnot(is.factor(x))
  if (length(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  x_split <- stri_split_regex(x, "\\.", simplify = TRUE)
  ## Which  column varies (recalling that after previous splits some columns may have no variance).
  ## try dataPreparation::whichAreConstant() or grab that function's C code in the future since it is much faster than below
  which_varies <- apply(x_split, 2, function(x) {
    length(unique(x))
  })
  ## Split on the first column with variance from the beginning
  split_on <- which(which_varies > 1)[1]
  if (is.na(split_on)) {
    group <- factor(rep(0, nrow(x)))
  } else {
    ## as.numeric of factors creates a vector that starts a 1
    group <- factor(as.numeric(x_split[, split_on] == x_split[1, split_on]))
  }
  return(group)
}

##' A set of pre-specified splits
##'
##' @param x is a data.table object where each column from 1 to k further divides the blocks.
##' Column one should be highest level and each other column should be nested
##' @param bid is not used
##' @export
splitSpecified <- function(bid, x) {
  stopifnot(is.data.table(x))
  if (nrow(x) == 2) {
    group <- factor(c(0, 1))
    return(group)
  }
  ## Which  column varies (recalling that after previous splits some columns may have no variance).
  which_varies <- x[, lapply(.SD, function(x) {
    length(unique(x))
  }), .SDcols = names(x), drop = TRUE]
  ## Split on the first column with variance from the beginning
  split_on <- names(which_varies)[as.vector(which_varies > 1)][1]
  if (is.na(split_on)) {
    group <- factor(rep(0, nrow(x)))
  } else {
    ## as.numeric of factors creates a vector that starts a 1
    group <- factor(as.numeric(x[, get(split_on)]) - 1)
  }
  return(group)
}
