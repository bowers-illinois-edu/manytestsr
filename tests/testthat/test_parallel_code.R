## Test and develop functions to adapt alpha levels as the tree grows
context("Comparing implementations of the basic distance based creation function for use in pIndepDist")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- TRUE
if (interactive) {
  library(here)
  library(data.table)
  source(here("tests/testthat", "make_test_data.R"))
  library(devtools)
  load_all() ## use  this during debugging
}


setDTthreads(1)
options(digits = 4)
## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)


test_length_fn <- function(obj1, obj2) {
  expect_equal(length(obj1), length(obj2))
}

test_contents_fn <- function(obj1, obj2) {
  expect_equal(obj1[[1]], obj2[[1]])
  expect_equal(obj1[[2]], obj2[[2]])
  expect_equal(obj1[[3]], obj2[[3]])
  expect_equal(obj1[[4]], obj2[[4]])
  expect_equal(obj1[[5]], obj2[[5]])
  expect_equal(obj1[[6]], obj2[[6]])
  # expect_equal(obj1[[7]], obj2[[7]])
  # expect_equal(obj1[[8]], obj2[[8]])
}

numcores <- parallel::detectCores(logical = FALSE)
numcores <- 4 # floor(cores/2)

test_that("All algorithmns give the same answers", {
  y <- rnorm(10)
  tmp1 <- dists_and_trans(y)
  tmp2 <- fast_dists_and_trans(y, Z = 1)
  #     List res = List::create(mndx, mndrx, maddx, maddrx, maxdx, maxdrx, zx, rankx);
  tmp3 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y, Z = 1)
  tmp5 <- fast_dists_by_unit_arma2_par(y, Z = 1, threads = numcores)
  #
  test_length_fn(tmp1, tmp2)
  test_length_fn(tmp1, tmp3)
  test_length_fn(tmp1, tmp4)
  test_length_fn(tmp1, tmp5)
  #
  test_contents_fn(tmp1, tmp2)
  test_contents_fn(tmp1, tmp3)
  test_contents_fn(tmp1, tmp4)
  test_contents_fn(tmp1, tmp5)
  #
  y <- rchisq(1000, df = 1)
  tmp1 <- dists_and_trans(y)
  tmp2 <- fast_dists_and_trans(y, Z = 1)
  tmp3 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y, Z = 1)
  tmp5 <- fast_dists_by_unit_arma2_par(y, Z = 1, threads = numcores)
  #
  test_length_fn(tmp1, tmp2)
  test_length_fn(tmp1, tmp3)
  test_length_fn(tmp1, tmp4)
  test_length_fn(tmp1, tmp5)
  #
  test_contents_fn(tmp1, tmp2)
  test_contents_fn(tmp1, tmp3)
  test_contents_fn(tmp1, tmp4)
  test_contents_fn(tmp1, tmp5)
})

## ## ### Check timing and core usage
## cores <- parallel::detectCores(logical = FALSE)
## set.seed(12345)
## y <- rchisq(5000, df = 1)
##
## res1 <- bench::mark(
##   main = dists_and_trans(y),
##   fastdandtrans = fast_dists_and_trans(y, Z = 1),
##   byunit_arma = fast_dists_and_trans_by_unit_arma(y, Z = 1),
##   byunit_arma2 = fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##   byunit_par = fast_dists_by_unit_arma2_par(y, Z = 1, threads = 7),
##   min_iterations = 100, max_iterations = 200, check = FALSE, filter_gc = FALSE
## )
## Seems like dists_and_trans is much faster with n=5000 (and also with n=100, n=1000, etc..)
# But it is using a ridiculous amount of memory (382MB versus 354K).
## The by unit approaches are so so much better given memory (for n=10000, like 1+GB versus 700 K)
# res1
# >  res1 (n=1000)
# # A tibble: 5 × 13
#      expression       min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory     time
#     <bch:expr>   <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>     <list>
#   1 main          21.2ms  21.8ms     45.3     15.4MB   0.453    100     1       2.2s <NULL> <Rprofmem> <bench_tm>
#   2 fastdandtra… 118.5ms   124ms      7.94    92.1MB   0.715    100     9      12.6s <NULL> <Rprofmem> <bench_tm>
#   3 byunit_arma  167.9ms 175.5ms      5.70    73.2KB   0        100     0      17.5s <NULL> <Rprofmem> <bench_tm>
#   4 byunit_arma2 166.1ms 169.6ms      5.87    73.2KB   0.0587   100     1        17s <NULL> <Rprofmem> <bench_tm>
#   5 byunit_par   167.1ms 170.2ms      5.88    73.2KB   0        100     0        17s <NULL> <Rprofmem> <bench_tm>
#   # ℹ 1 more variable: gc <list>
# (n=5000)
### A tibble: 5 × 13
##  expression         min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
##  <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>
## 1 main            89.6ms   92.4ms     10.7     61.3MB   8.66     100    81
## 2 fastdandtrans  525.6ms  602.5ms      1.70   367.2MB   8.12     100   477
## 3 byunit_arma    701.6ms  706.1ms      1.42   143.5KB   0        100     0
## 4 byunit_arma2   674.6ms  679.8ms      1.47   143.5KB   0.0147   100     1
## 5 byunit_par     673.6ms  678.6ms      1.47   143.5KB   0        100     0
##
## bench1 <- microbenchmark::microbenchmark(
##   main = dists_and_trans(y),
##   fastdandtrans = fast_dists_and_trans(y, Z = 1),
##   byunit_arma = fast_dists_and_trans_by_unit_arma(y, Z = 1),
##   byunit_arma2 = fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##   byunit_par = fast_dists_by_unit_arma2_par(y, Z = 1, threads = 7),
##   times = 100
## )
## bench1

## ## Look at how this works with the dpp data (moving most of this work to the test_pval_fns.R. Leaving this for now
## load(file = here::here("tests", "dpp_dat.rda"))
## table(dpp_dat$R01TMCRET, exclude = c())
## table(dpp_dat$R02TMCRET, exclude = c())
## ## Ranks more powerful. R02 easier to detect effects.
## pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pIndepDist(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
## pWilcox(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pWilcox(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
## pOneway(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pOneway(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)

## dpp_dat[,tmpy:=rchisq(.N,df=1)+trt*rnorm(.N),by=blockF]
## dpp_dat[,tmpy:=rnorm(.N)+trt*rnorm(.N),by=blockF]
## dpp_dat[, tmpy := R01TMCRET]
## tmpdat <- copy(dpp_dat)
## tmpdat[, dists_and_trans(tmpy), by = blockF]
## tmpdat[, c("mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY") := dists_and_trans(tmpy), by = blockF]
##
## tsmat <- as.matrix(tmpdat[, c("tmpy", "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")])
## cor(tsmat)
## solve(cor(tsmat))
##
## outcomenms <- c("tmpy", "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")
##
## test_fn <- function(indx) {
##   fmlatxt <- paste(paste(outcomenms[indx], collapse = "+"), "~trtF|blockF", collapse = "")
##   fmla <- as.formula(fmlatxt)
##   t1 <- independence_test(fmla, data = tmpdat, teststat = "maximum", distribution = asymptotic())
##   pvalue(t1)[[1]]
## }
##
## test_fn(indx = 1:length(outcomenms))
## ## Very little penalty for including highly correlated variables. (here even repeating a column)
## test_fn(indx = c(1, 1, 1:length(outcomenms)))
##
## res1 <- sapply(1:length(outcomenms), function(i) {
##   test_fn(indx = i)
## })
## names(res1) <- outcomenms
## res1
## ## Most power here is meandist, rank, std, raw --- mad is terrible as are the maximum distances.
## ##         tmpy       mndist  mndistRank0      maddist maddistRank0      maxdist maxdistRank0      zscoreY        rankY
## ##      0.05972      0.02646      0.01075      0.94665      0.96880      0.92711      0.70799      0.03700      0.03663
## ## Look at combinations of them
## all_but_one <- combn(9, 8)
## res2 <- apply(all_but_one, 2, function(i) {
##   test_fn(indx = i)
## })
## names(res2) <- apply(all_but_one, 2, function(idx) {
##   setdiff(1:ncol(all_but_one), idx)
## })
## res2
## ## Works well without the mad and max
## ## test_fn(indx=c(1,2,3,8,9))
