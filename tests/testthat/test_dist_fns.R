# Test and develop functions to produce distances between treated and control units

## The next lines are for use when creating the tests.
## Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(testthat)
  testthat::local_edition(3)
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
  testdist <- Rfast::Dist(xblah)
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
  manytestsr:::vecdist2(xblah)
}
vecdist3_fn <- function() {
  manytestsr:::vecdist3_arma(xblah)
}
vecdist_arma_fn <- function() {
  manytestsr:::vecdist_arma(xblah)
}
vecdist3_arma_fn <- function() {
  manytestsr:::vecdist3_arma(xblah)
}
## vecdist_rcpp_fn <- function() { vecdist_rcpp(xblah) }
## vecdist2_rcpp_fn <- function() { vecdist2_rcpp(xblah) }
basedist_fn <- function() {
  mat <- as.matrix(dist(xblah, diag = TRUE, upper = TRUE))
  dimnames(mat) <- NULL
  return(mat)
}
manhattan_fn <- function() {
  Rfast::Dist(xblah, method = "manhattan")
}
## manhattan2_fn <- function() {
##   manhattan_dist(xblah)
## }
## euc_dist_arma1_fn <- function() {
##   euc_dist_arma1(xblah)
## }

test_that("All distance matrix creation functions produce the same results", {
  expect_equal(byhand_fn(), Dist_fn())
  expect_equal(byhand_fn(), vecdist_fn())
  ## expect_equal(byhand_fn(), vecdist_rcpp_fn())
  ## expect_equal(byhand_fn(), vecdist2_rcpp_fn())
  expect_equal(byhand_fn(), vecdist2_fn())
  expect_equal(byhand_fn(), vecdist3_fn())
  expect_equal(byhand_fn(), vecdist3_arma_fn())
  expect_equal(byhand_fn(), vecdist_arma_fn())
  expect_equal(byhand_fn(), basedist_fn())
  expect_equal(byhand_fn(), manhattan_fn())
  # expect_equal(byhand_fn(), manhattan2_fn())
  # expect_equal(byhand_fn(), euc_dist_arma1_fn())
})

## Now test the score functions since we basically stopped making full
## euclidean distance matrices and have found it faster (and better for memory)
## to just work unit by unit in a loop in C

dists_and_trans_slow <- function(x) {
  n <- length(x)
  dx <- as.matrix(dist(x))
  dimnames(dx) <- c(NULL, NULL)
  rankx <- rank(x)
  dxRank0 <- as.matrix(dist(rankx)) # distance among the ranks
  dimnames(dxRank0) <- c(NULL, NULL)
  res <- list(
    mean_dist = colSums(dx) / (n - 1),
    mean_rank_dist = colSums(dxRank0) / (n - 1),
    max_dist = apply(dx, 2, max),
    rankY = rankx,
    tanhY = tanh(x)
  )
  return(res)
}

slow_r_fn <- function(x = xblah) {
  dists_and_trans_slow(x)
}

by_unit_fn <- function(x = xblah) {
  manytestsr:::fast_dists_and_trans_new(x)
}

r_fast_fn <- function(x = xblah) {
  dists_and_trans(x)
}

by_unit_fn_threads <- function(x = xblah) {
  manytestsr:::fast_dists_and_trans_new_omp(x, threads = 2)
}

just_four <- function(x = xblah) {
  manytestsr:::fast_dists_and_trans_nomax_hybrid(x)
}

test_that("All score creation functions produce same results", {
  expect_equal(by_unit_fn(), by_unit_fn_threads())
  expect_equal(slow_r_fn()[c("mean_dist", "mean_rank_dist", "rankY", "tanhY")], just_four())
  expect_equal(slow_r_fn(), r_fast_fn())
  expect_equal(by_unit_fn(), slow_r_fn())
  expect_equal(by_unit_fn(), r_fast_fn())
  expect_equal(manytestsr:::fast_dists_and_trans_new(xblah), fast_dists_and_trans_hybrid(xblah))
  expect_equal(manytestsr:::fast_dists_and_trans_new(xblah)[c("mean_dist", "mean_rank_dist", "rankY", "tanhY")], manytestsr:::fast_dists_and_trans_nomax_hybrid(xblah))
})


#### Benchmarks: We use the _hybrid functions in the pvalue code
##### Results from Macbook Pro M3 Max, 32GB RAM
# set.seed(12345)
# bench_n100 <- bench::mark(
#  cpp_new= fast_dists_and_trans_new(xblah),
#  cpp_hybrid= fast_dists_and_trans_hybrid(xblah),
#  cpp_hybrid_nomax= fast_dists_and_trans_nomax_hybrid(xblah),
#  r_fast= dists_and_trans(xblah),
#  r_slow= dists_and_trans_slow(xblah),
#  min_iterations = 100, max_iterations = 1000, check = FALSE, filter_gc = TRUE
# )
### Results from Macbook Pro M3 Max, 32GB RAM
# bench_n100 %>% arrange(median) %>% select(-result,-memory,-time,-gc)
## # A tibble: 5 × 9
##   expression            min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
##   <bch:expr>       <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
## 1 cpp_hybrid_nomax   12.8µs   14.2µs    69898.    3.31KB      0    1000     0     14.3ms
## 2 cpp_hybrid         18.4µs   20.2µs    49257.    4.14KB      0    1000     0     20.3ms
## 3 r_fast             21.8µs   27.8µs    35089.  161.31KB    247.    993     7     28.3ms
## 4 cpp_new            55.5µs   60.6µs    16326.    4.14KB      0    1000     0     61.3ms
## 5 r_slow            181.8µs  212.8µs     4691.  529.64KB     95.7   980    20    208.9ms
#
###
# set.seed(123455)
# n <- 1000
# xblah_n1000 <- rchisq(n, df = 1) + rnorm(n)
# set.seed(12345)
# bench_n1000 <- bench::mark(
#  cpp_new= fast_dists_and_trans_new(xblah_n1000),
#  cpp_hybrid= fast_dists_and_trans_hybrid(xblah_n1000),
#  cpp_hybrid_nomax= fast_dists_and_trans_nomax_hybrid(xblah_n1000),
#  r_fast= dists_and_trans(xblah_n1000),
#  r_slow= dists_and_trans_slow(xblah_n1000),
#  min_iterations = 100, max_iterations = 1000, check = FALSE, filter_gc = TRUE
# )
# bench_n1000 %>% arrange (median) %>% select(-result,-memory,-time,-gc)
## # A tibble: 5 × 9
##   expression            min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
##   <bch:expr>       <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
## 1 cpp_hybrid_nomax 118.53µs 130.38µs    7669.     31.4KB     0     1000     0    130.4ms
## 2 cpp_hybrid       173.27µs 191.76µs    5223.     39.3KB     5.23   999     1    191.3ms
## 3 r_fast             1.62ms   1.91ms     516.     15.3MB   297.      99    57    191.7ms
## 4 cpp_new            4.42ms   4.56ms     219.     39.3KB     2.02   108     1    494.1ms
## 5 r_slow            13.74ms  13.74ms      72.8    49.8MB  8659.       1   119     13.7ms
#
## set.seed(123455)
## n <- 10000
## xblah_n10000 <- rchisq(n, df = 1) + rnorm(n)
## set.seed(12345)
## bench_n10000 <- bench::mark(
##  cpp_sort= fast_dists_and_trans_sort(xblah_n10000),
##  cpp_new= fast_dists_and_trans_new(xblah_n10000),
##  cpp_hybrid= fast_dists_and_trans_hybrid(xblah_n10000),
##  r_fast= dists_and_trans(xblah_n10000),
##  r_slow= dists_and_trans_slow(xblah_n10000),
##  min_iterations = 100, max_iterations = 1000, check = TRUE, filter_gc = TRUE
## )
#### Results from Macbook Pro M3 Max, 32GB RAM
## bench_n10000 %>% arrange (median) %>% select(-result,-memory,-time,-gc)
### # A tibble: 5 × 9
###   expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time
###   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm>
### 1 cpp_hybrid   1.93ms   2.06ms   477.     390.86KB     0      239     0   500.72ms
### 2 r_fast     265.15ms  283.4ms     3.44     1.49GB     3.68   100   107     29.05s
### 3 cpp_sort   398.79ms 401.83ms     2.47   390.86KB     0      100     0     40.56s
### 4 cpp_new    466.79ms 469.19ms     2.11   390.86KB     0      100     0     47.32s
### 5 r_slow        1.73s    1.79s     0.557    4.84GB     2.38   100   428      2.99m


#

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
