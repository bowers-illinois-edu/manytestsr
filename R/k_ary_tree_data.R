#' Example tree-testing data
#'
#' Synthetic data used in examples and tests for \pkg{manytestsr}.
#'
#' Generated once and stored as package data to avoid slow/regenerative
#' simulation in examples/tests. The full recipe lives in
#' \file{data-raw/make_k_ary_tree_data.R}.
#'
#' @docType data
#' @name idt
#' @aliases bdt1
#' @format
#' \describe{
#'   \item{idt}{Individual-level data.table with outcomes and assignments.}
#'   \item{bdt1}{Block-level data.table containing summaries per block.}
#' }
#' @source See \file{data-raw/make_k_ary_tree_data.R}; produced with
#'   \pkg{TreeTestSim} (≥ X.Y.Z) and a fixed RNG seed on 2025-11-08.
#' @references Bowers et al., \pkg{TreeTestSim}, \pkg{manytestsr}.
#' @keywords datasets
#' @examples
#' str(idt)
#' str(bdt1)
"idt"

#' @rdname idt
"bdt1"
