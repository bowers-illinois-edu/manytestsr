## Test Statistics

#' The combined distances test statistic

#' @examples
#' \dontrun{
#' # Requires dat object to be defined in environment
#' tstat_dt <- combined_distances_tstat()
#'
#' splitCluster("blockF", x)
#' ## clus <- tryCatch(KMeans_rcpp(as.matrix(x), clusters = 2, num_init = 2,
#' ##        initializer = "optimal_init")$clusters, error = function(e) {
#' ##    kmeans(x, centers = 2)$cluster})
#'
#' ## Approach 2:
#' clus <- Ckmeans.1d.dp(x, k = 2)$cluster
#' group <- factor(as.numeric(clus == 1))
#' }

#' @export
combined_distances_tstat <- function(fmla = Y ~ trtF | blockF, distfn = fast_dists_and_trans_hybrid) {
  force(distfn)
  stopifnot(inherits(dat, "data.table"))
  fmla_vars <- all.vars(fmla)
  theresponse <- fmla_vars[attr(terms(fmla), "response")]
  thetreat <- fmla_vars[[2]]

  thedat <- copy(dat)

  ### These must match the names of the functions used in src/dists_and_trans.cpp
  outcome_names <- c(theresponse, "mean_dist", "mean_rank_dist", "max_dist", "rankY", "tanhY")

  if (length(fmla_vars) == 3) {
    theblock <- fmla_vars[[3]]
    thedat[, outcome_names[-1] := distfn(get(theresponse)), by = get(theblock)]
  } else {
    theblock <- NULL
    # This next is faster than doing it in two lines
    thedat[, outcome_names[-1] := distfn(get(theresponse))]
  }

  # If one of the test statistics is constant, drop it.
  # https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
  anyconstant_cols <- dataPreparation::which_are_constant(thedat[, .SD, .SDcols = outcome_names], verbose = FALSE)
  if (length(anyconstant_cols) > 0) {
    outcome_names <- outcome_names[-anyconstant_cols]
  }
  #    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, "|", theblock, sep = "")
  #    newfmla_text <- paste(paste(outcome_names, collapse = "+"), "~", thetreat, sep = "")
  #  newfmla <- as.formula(newfmla_text)

  if (!is.null(theblock)) {
    cols_to_return <- c(outcome_names, theresponse, thetreat, theblock)
  } else {
    cols_to_return <- c(outcome_names, theresponse, thetreatD)
  }
  return(thedat[, .SD, .SDcols = cols_to_return])
}
