## Alpha adjustment functions for sequential testing based FDR control
## For now using the onlineFDR package functions


##' Alpha adjustment function: Alpha Investing
##'
##' These function accept a vector of p-values and vector indicating batch of p-values
##' It produces a vector of alpha values.
##'
##' @param pval Is a numeric vector of p-values
##' @param batch A character or factor or even numeric vector indicating grouping among the p-values
##' and the order of the p-values. The first element pval and first element of batch should be the first p-value
##' produced in the tree and the subsequent p-values should be produced on the children of that root (or children of the children of the root, etc.)
##' @param nodesize Is vector indicating the information used to create each p-value. For example, it could be the number of observations, or a weighted version.
##' @return A vector of alpha values
##' @importFrom onlineFDR Alpha_investing
##' @import data.table
##' @export
alpha_investing <- function(pval, batch, nodesize) {
  stopifnot(length(pval) == length(nodesize))
  stopifnot(length(batch) == length(nodesize))
  ## Later just use batch when we fork onlineFDR.
  ## For now the onlineFDR functions want dates to indicate batches. A hassle
  ### Since we sometimes have more than 31 batches we just convert into days since 01-01-0001.
  if (is.numeric(batch)) {
    batchuniq <- as.Date(batch, origin = "0001-01-01")
  } else {
    batchuniq <- as.Date(as.numeric(as.factor(batch)), origin = "0001-01-01")
  }
  ## id is required by Alpha_investing although it is meaningless for us
  thedf <- data.table(id = seq(1, length(pval)), pval = pval, date = batchuniq, nodesize = nodesize)
  ## Order tests by nodesize within batch: bigger nodes come first, they should have more power and thus smaller alpha
  thedf <- thedf[order(date, -nodesize), ]
  res <- Alpha_investing(d = thedf, random = FALSE)
  return(res$alphai[thedf$id])
}

##' Alpha adjustment function: SAFFRON
##'
##' These function accept a vector of p-values and vector indicating batch of p-values
##' It produces a vector of alpha values.
##'
##' @param pval Is a numeric vector of p-values
##' @param batch A character or factor or even numeric vector indicating grouping among the p-values
##' and the order of the p-values. The first element pval and first element of batch should be the first p-value
##' produced in the tree and the subsequent p-values should be produced on the children of that root (or children of the children of the root, etc.)
##' @param nodesize Is vector indicating the information used to create each p-value. For example, it could be the number of observations, or a weighted version.
##' @return A vector of alpha values
##' @importFrom onlineFDR SAFFRON
##' @export
alpha_saffron <- function(pval, batch, nodesize) {
  stopifnot(length(pval) == length(nodesize))
  stopifnot(length(batch) == length(nodesize))
  ## Later just use batch. For now the onlineFDR functions want dates to indicate batches
  if (is.numeric(batch)) {
    batchuniq <- as.Date(batch, origin = "0001-01-01")
  } else {
    batchuniq <- as.Date(as.numeric(as.factor(batch)), origin = "0001-01-01")
  }
  ## id is required by Alpha_investing although it is meaningless for us
  thedf <- data.table(id = seq(1, length(pval)), pval = pval, date = batchuniq, nodesize = nodesize)
  ## Order tests by nodesize within batch: bigger nodes come first, they should have more power and thus smaller alpha
  thedf <- thedf[order(date, -nodesize), ]
  res <- SAFFRON(d = thedf, random = FALSE)
  return(res$alphai[thedf$id])
}
