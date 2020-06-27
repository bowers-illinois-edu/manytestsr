# Distance Creation Functions

#' Transformation Helper Function: Calculate Interunit  Distances
#'
#' These functions a  vector and produce a data.frame containing  distances.
#'
#' @param x A  numeric vector
#' @return A numeric vector of the sums  of  centered distances from  each unit the other units
#' @export
distsums_raw <- function(x) {
  # An  improvisation off the work of Rizzo  and her collaborators on e-statistics
  # Using their Akl function in dcov.c  in the energy package
  dx <- as.matrix(dist(x))
  # dx_centered <- dx - colMeans(dx)  -  rowMeans(dx) + mean(dx)
  # res <- data.frame(Ytrans =.rowSums(dx_centered,m=nrow(dx),n=ncol(dx)))
  res <- data.frame(Ytrans = .rowSums(dx, m = nrow(dx), n = ncol(dx)))
  return(res)
}

#' @export
distsums_abs <- function(x) {
  dx <- as.matrix(dist(x))
  dx_centered <- dx - colMeans(dx) - rowMeans(dx) + mean(dx)
  res <- data.frame(Ytrans = .rowSums(abs(dx_centered), m = nrow(dx), n = ncol(dx)))
  return(res)
}

#' @export
distsums_sq <- function(x) {
  dx <- as.matrix(dist(x))
  res <- data.frame(Ytrans = .rowSums(dx^2, m = nrow(dx), n = ncol(dx)))
  return(res)
}

#' @export
distsums_stf <- function(x) {
  dx <- as.matrix(dist(x))
  n <- nrow(dx)
  dx_se <- sqrt(apply(dx, 1, var) / n)
  res <- data.frame(Ytrans = .rowMeans(dx, m = nrow(dx), n = ncol(dx)) / dx_se)
  return(res)
}


#' Outcome e-distances between treatment arms
#'
#' @param x Is a numeric vector (the  outcome variable)
#' @param Z is a binary numeric vector or a  factor vector with  only  two values
#' @return Vector with  the individual level components of the energy distances between each unit and the units in  the control condition
#' @importFrom Rfast Dist colsums rowsums Order
#' @export
edisti <- function(x, Z) {
  dx <- vecdist(x) ## could be dx^2 for anova style
  if (is.factor(Z)) {
    Z <- as.numeric(levels(Z))[Z] - 1
    stopifnot(all(Z %in% c(0, 1)))
  }
  nT <- sum(Z)
  nC <- sum(1 - Z)
  ix <- 1:length(x)
  names(x) <- ix
  names(Z) <- ix
  dimnames(dx) <- list(ix, ix)
  Tix <- ix[Z == 1]
  Cix <- ix[Z == 0]
  dxT <- dx[Tix, Tix]
  dxC <- dx[Cix, Cix]
  dxTC <- dx[Tix, Cix]
  thew <- (nT * nC) / (nT + nC)
  # TCi <- Rfast::colsums(dxTC)
  TCiC <- colsums(dxTC) # results in order of Cix
  names(TCiC) <- Cix # This just for checking ordering. Could do it more efficiently with stats::colSums()
  TCiT <- rowsums(dxTC) # results in order of Tix
  names(TCiT) <- Tix
  Ti <- colsums(dxT) # order of Tix, the dxT and dxC mats will be square so col or row is same
  names(Ti) <- Tix
  Ci <- colsums(dxC) # order of Cix, the dxT and dxC mats will be square so col or row is same
  names(Ci) <- Cix

  # This next matches edist()
  # e_mine2 <- thew * ( ( sum(TCiC)/(nT*nC) - sum(Ci)/(nC*nC) ) + ( sum(TCiT)/(nT*nC) - sum(Ti)/(nT*nT) ) )

  eCi <- ((TCiC / (nT * nC)) - (Ci / (nC * nC))) # assuming same order of TCiC and Ci
  eTi <- ((TCiT / (nT * nC)) - (Ti / (nT * nT))) # assuming same order of TCiT and Ti

  ei <- (c(eCi, eTi) * thew)[Order(c(Cix, Tix))]
  stopifnot(all.equal(names(ei), names(Z)))
  stopifnot(all.equal(names(ei), names(x)))
  # These also work
  # e_mine3 <- thew*sum(ei)
  # e_mine4 <- sum(ei* thew )
  # So e per person is something like this
  # ei2 <- c(eCi*thew,eTi*thew)
  # This also works:
  # e_mine5 <- sum(ei2)
  # m11 <- sum(dst[ii, ii]) / (n1 * n1)
  # m22 <- sum(dst[jj, jj]) / (n2 * n2)
  # m12 <- sum(dst[ii, jj]) / (n1 * n2)
  # e[i, j] <- e[j, i] <- w * ((m12 + m12) - (m11 + m22))
  # res <- w * ( (c_t_dists - c_c_dists) + (c_t_dists - t_t_dists) )
  return(ei)
}

#' Outcome distances between treatment arms and transformations of distances
#'
#' @param x Is a numeric vector (the  outcome variable)
#' @param Z is just a placeholder and not used but is a part of the api for distance functions
#' @return A list  for inclusion in a data.table with  distances between each unit and other units as well as some transformations of those distances
#' @importFrom Rfast Dist rowmeans rowMads rowMaxs mahala cova
#' @export
dists_and_trans <- function(x, Z) {
  # I bet we can speed this up by just doing it in cpp
  ## if(!is.numeric(x){x <- as.numeric(x)}
  dx <- vecdist(x)
  dxRank0 <- vecdist(Rank(x)) # distance among the ranks
  mnx <- fastMean(x)
  xmat <- matrix(x, ncol = 1)
  sigx <- cova(xmat)
  res <- list(
    mndist = rowmeans(dx),
    mndistRank0 = rowmeans(dxRank0),
    maddist = rowMads(dx),
    maddistRank0 = rowMads(dxRank0),
    maxdist = rowMaxs(dx, value = TRUE),
    maxdistRank0 = rowMaxs(dxRank0, value = TRUE),
    mhdist = zscore_vec(x)
  )
  return(res)
}


#
# #' OLD: Outcome distances between treatment arms
# #'
# #' @param x Is a numeric vector (the  outcome variable)
# #' @param Z is a binary numeric vector or a  factor vector with  only  two values
# #' @return A list  for inclusion in a data.table with  distances between each unit and other units
# olddistfns <- function(x, Z) {
#   # allow distance among the ranks?
#   # dx <- as.matrix(stats::dist(x)) ## could be dx^2 for anova style
#   dx <- Rfast::Dist(x)
#   if (is.factor(Z)) {
#     Z <- as.numeric(levels(Z))[Z] - 1
#     stopifnot(all(Z %in% c(0, 1)))
#   }
#   ntncmat <- diag(1 / c(sum(Z), sum(1 - Z)))
#   M <- cbind(other = (1 - Z), own = Z)
#   # sums of distances from each to controls (column 1) and to treated (column 2)
#   # ctrl_trt_sums <- t(M) %*% dx
#   # ctrl_trt_sums <- dx %*% M # same
#   ctrl_trt_sums <- eigenMapMatMult(dx, M)
#   # ctrl_trt_sums has two columns. First column is sum of distances to controls. Second column is sum of distances to treated.
#   # This next is not the same as energy distance but attempts to make something that could add up to it or
#   # add up to something proportionate to it
#   ctrl_trt_avg_sums <- ctrl_trt_sums %*% ntncmat
#   cross_avg_dists <- (ctrl_trt_avg_sums * (1 - M)) %*% c(1, 1) # n x 2 * n x 2 %*% 2 x 1
#   within_avg_dists <- (ctrl_trt_avg_sums * M) %*% c(1, 1)
#   indiv_diff <- cross_avg_dists - within_avg_dists
#   indiv_rat <- cross_avg_dists / within_avg_dists
#   # Now with ranks
#   dxRank0 <- Rfast::Dist(Rank(x)) # distance among the ranks
#   ctrl_trt_sums_Rank0 <- eigenMapMatMult(dxRank0, M)
#   ctrl_trt_avg_sums_Rank0 <- ctrl_trt_sums_Rank0 %*% ntncmat
#   cross_avg_dists_Rank0 <- (ctrl_trt_avg_sums_Rank0 * (1 - M)) %*% c(1, 1) # n x 2 * n x 2 %*% 2 x 1
#   within_avg_dists_Rank0 <- (ctrl_trt_avg_sums_Rank0 * M) %*% c(1, 1)
#   indiv_diff_Rank0 <- cross_avg_dists_Rank0 - within_avg_dists_Rank0
#   indiv_rat_Rank0 <- cross_avg_dists_Rank0 / within_avg_dists_Rank0
#   dxRank1 <- Rfast::rowRanks(dx) # ranks of the distances
#   diag(dxRank1) <- 0
#   ctrl_trt_sums_Rank1 <- eigenMapMatMult(dxRank1, M)
#   ctrl_trt_avg_sums_Rank1 <- ctrl_trt_sums_Rank1 %*% ntncmat
#   cross_avg_dists_Rank1 <- (ctrl_trt_avg_sums_Rank1 * (1 - M)) %*% c(1, 1) # n x 2 * n x 2 %*% 2 x 1
#   within_avg_dists_Rank1 <- (ctrl_trt_avg_sums_Rank1 * M) %*% c(1, 1)
#   indiv_diff_Rank1 <- cross_avg_dists_Rank1 - within_avg_dists_Rank1
#   indiv_rat_Rank1 <- cross_avg_dists_Rank1 / within_avg_dists_Rank1
#   # Returns a list because of use within data.table
#   res <- list(
#     # indiv_e_scordiff = drop ( indiv_diff ),
#     mndist = Rfast::rowmeans(dx), # trying something even simpler
#     mndistRank0 = Rfast::rowmeans(dxRank0),
#     # mndistRank1 = Rfast::rowmeans(dxRank1),
#     maddist = Rfast::rowMads(dx),
#     maddistRank0 = Rfast::rowMads(dxRank0),
#     # maddistRank1 = Rfast::rowMads(dxRank1),
#     maxdist = Rfast::rowMaxs(dx, value = TRUE), # trying something even simpler
#     maxdistRank0 = Rfast::rowMaxs(dxRank0, value = TRUE) # ,
#     # maxdistRank1 = Rfast::rowMaxs(dxRank1,value=TRUE)
#   )
#   return(res)
# }
