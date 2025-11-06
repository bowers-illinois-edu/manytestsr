# Distance Creation Functions


#' Outcome distances and transformations
#'
#' @param x Is a numeric vector (the  outcome variable)
#' @return A list  for inclusion in a data.table with  distances between each unit and other units as well as some transformations of those distances
#' @examples
#' # Example with continuous outcome data
#' outcome <- c(2.1, 5.3, 1.8, 7.2, 3.4, 4.6)
#' 
#' # Compute distance-based transformations
#' dist_results <- dists_and_trans(outcome)
#' 
#' # View components
#' str(dist_results)
#' print(dist_results$mean_dist)     # Mean distance to all other units
#' print(dist_results$rankY)         # Ranks of original values
#' print(dist_results$tanhY)         # Hyperbolic tangent transformation
#' 
#' @importFrom Rfast Dist rowmeans rowMads rowMaxs Rank vecdist rowmeans
#' @export
dists_and_trans <- function(x) {
  n <- length(x)
  dx <- Rfast::vecdist(x)
  rankx <- Rfast::Rank(x)
  dxRank0 <- Rfast::vecdist(rankx) # distance among the ranks
  res <- list(
    ## Don't include the diagonals in the means
    mean_dist= Rfast::colsums(dx)/(n-1),
    mean_rank_dist= Rfast::colsums(dxRank0)/(n-1),
    max_dist = Rfast::colMaxs(dx, value = TRUE),
    rankY = rankx,
    tanhY = tanh(x)
  )
  return(res)
}

#' Outcome distances and transformations: C++ OpenMP Parallel version
#'
#' @param threads Is an integer with the number of cores to use.
#' @return A distance creation function taking x (a numeric vector, usually the outcome) and a variable Z which is not used.
#' @export
fast_dists_and_trans_new_parallel <- function(threads) {
  force(threads)
  return(function(x) {
    fast_dists_and_trans_new_omp(x, threads = threads)
  })
}




#' Outcome e-distances between treatment arms
#'
#' @param x Is a numeric vector (the  outcome variable)
#' @param Z is a binary numeric vector or a  factor vector with  only  two values
#' @return Vector with  the individual level components of the energy distances between each unit and the units in  the control condition
#' @examples
#' # Example with treatment and control groups
#' outcome <- c(2.1, 5.3, 1.8, 7.2, 3.4, 4.6, 6.1, 2.8)
#' treatment <- c(0, 1, 0, 1, 0, 1, 1, 0)  # Binary treatment indicator
#'
#' # Compute energy distances
#' e_dists <- edisti(outcome, treatment)
#' print(e_dists)
#'
#' \dontrun{
#' # With factor treatment variable
#' treatment_factor <- factor(treatment, labels = c("Control", "Treatment"))
#' e_dists2 <- edisti(outcome, treatment_factor)
#' print(e_dists2)
#' }
#' 
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
