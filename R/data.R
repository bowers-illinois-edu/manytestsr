#' Example Data for Block-Randomized Experiment
#'
#' A simulated dataset demonstrating a block-randomized experiment with
#' two outcome variables. This dataset is used throughout the package examples
#' to illustrate the use of various testing and splitting procedures.
#'
#' @format A data.table with 1268 rows and 9 variables:
#' \describe{
#'   \item{id}{Integer. Unique identifier for each unit.}
#'   \item{year}{Integer. Year indicator (1 or 3).}
#'   \item{trt}{Integer. Binary treatment indicator (0 = control, 1 = treated).}
#'   \item{Y1}{Numeric. First outcome variable.}
#'   \item{Y2}{Numeric. Second outcome variable.}
#'   \item{trtF}{Factor. Treatment indicator as a factor with levels "0" and "1".}
#'   \item{place_year_block}{Character. Combined identifier for place, year, and block.}
#'   \item{place}{Character. Place/location identifier (A, B, C, etc.).}
#'   \item{blockF}{Factor. Block identifier with 44 levels (B080 through B123).}
#' }
#'
#' @source Simulated data for demonstration purposes
#'
#' @examples
#' data(example_dat, package = "manytestsr")
#' head(example_dat)
#' table(example_dat$trt)
#' table(example_dat$blockF)
"example_dat"
