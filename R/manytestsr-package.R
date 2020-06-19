#' @useDynLib manytestsr, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.formula dist kmeans median na.omit p.adjust rnorm runif sd terms var
NULL

# a fix for R check warnings from https://github.com/STAT545-UBC/Discussion/issues/451
# see for example https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
