# Distance Creation Functions


#' Outcome distances and transformations
#'
#' @param x Is a numeric vector (the  outcome variable)
#' @param Z is just a placeholder and not used but is a part of the api for distance functions
#' @return A list  for inclusion in a data.table with  distances between each unit and other units as well as some transformations of those distances
#' @importFrom Rfast Dist rowmeans rowMads rowMaxs Rank vecdist rowmeans
#' @export
dists_and_trans <- function(x, Z) {
  dx <- Rfast::vecdist(x)
  rankx <- Rfast::Rank(x)
  dxRank0 <- Rfast::vecdist(rankx) # distance among the ranks
  mhdist0 <- zscore_vec(x) ## calling it mhdist but really just zscore
  mhdist <- ifelse(is.na(mhdist0)|is.nan(mhdist0),0,mhdist0)
  res <- list(
    mndist = Rfast::rowmeans(dx),#as.numeric(fastrowMeans(dx)),
    mndistRank0 = Rfast::rowmeans(dxRank0),#as.numeric(fastrowMeans(dxRank0)),
   # maddist = Rfast::rowMads(dx),
   # maddistRank0 = Rfast::rowMads(dxRank0),
    maxdist = Rfast::rowMaxs(dx,value=TRUE), #as.numeric(fastrowMaxs2(dx)),
    maxdistRank0 = Rfast::rowMaxs(dxRank0,value=TRUE), #as.numeric(fastrowMaxs2(dxRank0)),
    mhdist = mhdist,
    rankx = rankx
  )
  return(res)
}

#' Outcome distances and transformations: C++ OpenMP Parallel version
#'
#' @param threads Is an integer with the number of cores to use.
#' @return A distance creation function taking x (a numeric vector, usually the outcome) and a variable Z which is not used.
#' @export
fast_dists_by_unit_arma_parR <- function(threads) {
    force(threads)
    return(function(x,Z){ fast_dists_by_unit_arma2_par(x,Z,threads=threads) })
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
