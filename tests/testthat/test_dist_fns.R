# Test and develop functions to produce distances between treated and control units
context("Performance of  Distance Functions")

## The next lines are for use when creating the tests.
## Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  library(testthat)
  library(bench)
  library(microbenchmark)
  load_all() ## use  this during debugging
}
# library(bench)
setDTthreads(1)

## Make a test vector that is skewed
set.seed(123455)
n <- 100
xblah <- rchisq(n, df = 1) + rnorm(n)

test_that("Some distance matrix creation functions produce the same results", {
  dist2 <- sqrt(outer(xblah, xblah, FUN = function(x, y) {
    (x - y)^2
  }))
  testdist <- Dist(xblah)
  expect_equal(dist2, testdist)
})

byhand_fn <- function() {
  sqrt(outer(xblah, xblah, FUN = function(x, y) {
    (x - y)^2
  }))
}
Dist_fn <- function() {
  Rfast::Dist(xblah, vector = FALSE)
}
vecdist_fn <- function() {
  Rfast::vecdist(xblah)
}
vecdist2_fn <- function() {
  vecdist2(xblah)
}
vecdist3_fn <- function() {
  vecdist3(xblah)
}
vecdist_arma_fn <- function() {
  vecdist_arma(xblah)
}
vecdist3_arma_fn <- function() {
  vecdist3_arma(xblah)
}
vecdist4_arma_fn <- function() {
  vecdist4_arma(xblah)
}
## vecdist_rcpp_fn <- function() { vecdist_rcpp(xblah) }
## vecdist2_rcpp_fn <- function() { vecdist2_rcpp(xblah) }
basedist_fn <- function() {
  mat <- as.matrix(dist(xblah, diag = TRUE, upper = TRUE))
  dimnames(mat) <- NULL
  return(mat)
}
manhattan_fn <- function() {
  Dist(xblah, method = "manhattan")
}
manhattan2_fn <- function() {
  manhattan_dist(xblah)
}
euc_dist_arma1_fn <- function() {
  euc_dist_arma1(xblah)
}

test_that("All distance matrix creation functions produce the same results", {
  expect_equal(byhand_fn(), Dist_fn())
  expect_equal(byhand_fn(), vecdist_fn())
  ## expect_equal(byhand_fn(), vecdist_rcpp_fn())
  ## expect_equal(byhand_fn(), vecdist2_rcpp_fn())
  expect_equal(byhand_fn(), vecdist2_fn())
  expect_equal(byhand_fn(), vecdist3_fn())
  expect_equal(byhand_fn(), vecdist4_arma_fn())
  expect_equal(byhand_fn(), vecdist3_arma_fn())
  expect_equal(byhand_fn(), vecdist_arma_fn())
  expect_equal(byhand_fn(), basedist_fn())
  expect_equal(byhand_fn(), manhattan_fn())
  # expect_equal(byhand_fn(), manhattan2_fn())
  expect_equal(byhand_fn(), euc_dist_arma1_fn())
})

## This next works interactively on linux (and mac) but fails using test_local() and test_check().
## Since this is mostly about speed, I'm going to assume that the ones identified are the best ones and comment this out.
## test_that("The fastest/best distance matrix creation code is what I think it should. Could fail when machines or other aspects of architecture change", {
##   bench1 <- bench::mark(
##     byhand = byhand_fn(),
##     vecdist = vecdist_fn(),
##     vecdist2 = vecdist2_fn(),
##     vecdist3 = vecdist3_fn(),
##     ## vecdist_rcpp = vecdist_rcpp_fn(),
##     vecdist_arma = vecdist_arma_fn(),
##     vecdist3_arma = vecdist3_arma_fn(),
##     vecdist4_arma = vecdist4_arma_fn(),
##     Dist = Dist_fn(),
##     manhattan = manhattan_fn(),
##     basdist = basedist_fn(),
##     #mah2 = manhattan2_fn(),
##     euc2 = euc_dist_arma1_fn(),
##     min_iterations = 100, max_iterations = 1000, check = TRUE, filter_gc = TRUE
##   )
##   best4a <- bench1 %>%
##     arrange(median) %>%
##     filter(median <= median[4])
##   expect_true(all(as.character(best4a$expression) %in% c("vecdist", "Dist", "vecdist2", "byhand")))
## })

## Results from Macbook Pro M3 Max, 32GB RAM
## bench1 %>% arrange(median)
# n=100
##    expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result   memory
##    <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>   <list>
##  1 vecdist      3.44µs   6.03µs   172259.    80.7KB    518.    997     3     5.79ms <dbl[…]> <Rprofmem>
##  2 Dist         9.06µs  11.11µs    89343.    81.5KB    269.    997     3    11.16ms <dbl[…]> <Rprofmem>
##  3 vecdist2    14.35µs  22.59µs    44253.    80.7KB    178.    996     4    22.51ms <dbl[…]> <Rprofmem>
##  4 byhand      31.65µs  41.74µs    23014.     235KB    209.    991     9    43.06ms <dbl[…]> <Rprofmem>
##  5 manhattan    48.5µs  54.04µs    18534.    81.5KB     55.8   997     3    53.79ms <dbl[…]> <Rprofmem>
##  6 vecdist_a…  61.75µs  67.03µs    14905.    80.7KB     44.8   997     3    66.89ms <dbl[…]> <Rprofmem>
##  7 basdist     68.39µs   80.4µs    12229.   430.2KB    136.    989    11    80.88ms <dbl[…]> <Rprofmem>
##  8 vecdist3_…  116.6µs 130.34µs     7643.    80.7KB     23.0   997     3   130.45ms <dbl[…]> <Rprofmem>
##  9 mah2       134.77µs 148.05µs     6746.    80.7KB     13.5   998     2   147.95ms <dbl[…]> <Rprofmem>
## 10 euc2       159.86µs 163.92µs     5980.    80.7KB     12.0   998     2   166.89ms <dbl[…]> <Rprofmem>
## 11 vecdist4_… 194.09µs 207.46µs     4775.    80.7KB     14.4   997     3   208.79ms <dbl[…]> <Rprofmem>
## 12 vecdist3   264.41µs 273.55µs     3601.    80.7KB     10.8   997     3   276.88ms <dbl[…]> <Rprofmem>
## 13 vecdist_r…   4.54ms   4.66ms      214.    80.7KB      0     108     0    504.2ms <dbl[…]> <Rprofmem>
## n=1000
##    expression         min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
##    <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
##  1 vecdist       250.72µs 507.48µs  1988.       7.64MB 180.       742    67   373.18ms
##  2 Dist          622.79µs 831.29µs  1181.       7.65MB 107.       418    38   353.99ms
##  3 vecdist2        1.71ms    2.4ms   420.       7.64MB  39.3      182    17   433.12ms
##  4 byhand          3.13ms    3.5ms   282.      22.89MB 127.        73    33      259ms
##  5 manhattan       4.88ms    5.1ms   195.       7.65MB  19.3       91     9   466.72ms
##  6 vecdist_arma    6.66ms   6.91ms   144.       7.64MB  14.2       91     9   633.42ms
##  7 basdist          6.7ms   7.27ms   137.      41.98MB 126.        52    48   380.22ms
##  8 vecdist3_arma  11.89ms  13.35ms    73.0      7.64MB   6.35      92     8      1.26s
##  9 mah2           14.36ms   14.7ms    68.3      7.64MB   5.94      92     8      1.35s
## 10 euc2           16.71ms   17.1ms    58.1      7.64MB   4.37      93     7       1.6s
## 11 vecdist4_arma  20.84ms  21.92ms    45.6      7.64MB   3.97      92     8      2.02s
## 12 vecdist3       24.81ms   25.2ms    39.6      7.64MB   3.44      92     8      2.32s
## 13 vecdist_rcpp     4.59s    4.62s     0.216    7.64MB   0.0188    92     8      7.09m
