#' Example tree-testing data
#'
#' Synthetic data used in examples and tests for \pkg{manytestsr}.
#'
#' Generated once and stored as package data to avoid slow/regenerative
#' simulation in examples/tests. The full recipe lives in
#' \file{data-raw/make_k_ary_tree_data.R}.
#'
#' @format
#' \describe{
#'   \item{idt}{individual level data.table}
#'   \item{bt1}{block level data.table}
#' }
#' @source See \file{data-raw/make_k_ary_tree_data.R}; produced with
#'   \pkg{TreeTestSim} (≥ X.Y.Z) and a fixed RNG seed on 2025-11-08.
#' @references Bowers et al., \pkg{TreeTestSim}, \pkg{manytestsr}.
#' @keywords datasets
#' @examples
#' str(idt)
#' str(bdt1)
"idt"
"bdt1"
