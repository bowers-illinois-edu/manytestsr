## Test and develop functions to adapt alpha levels as the tree grows
context("Comparing implementations of the basic distance based creation function for use in pIndepDist")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
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
  lenobj1 <- length(obj1)
  lenobj2 <- length(obj2)
  stopifnot(all.equal(lenobj1, lenobj2))
  for (i in 1:lenobj1) {
    expect_equal(obj1[[i]], obj2[[i]])
  }
}

numcores <- parallel::detectCores(logical = FALSE)
numcores <- 4 # floor(cores/2)

test_that("All algorithms give the same answers", {
  y <- rnorm(10)
  tmp1 <- dists_and_trans(y)
  tmp2 <- fast_dists_and_trans(y)
  tmp3 <- fast_dists_and_trans_by_unit_arma(y)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y)
  tmp5 <- fast_dists_and_trans_by_unit_arma2_par(y, threads = numcores)
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
  tmp2 <- fast_dists_and_trans(y)
  tmp3 <- fast_dists_and_trans_by_unit_arma(y)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y)
  tmp5 <- fast_dists_and_trans_by_unit_arma2_par(y, threads = numcores)
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

## ## ## ### Check timing and core usage
## cores <- parallel::detectCores(logical = FALSE)
## set.seed(12345)
## y <- rchisq(100, df = 1)
##
## res1 <- bench::mark(
##   main = dists_and_trans(y),
##   fastdandtrans = fast_dists_and_trans(y),
##   byunit_arma = fast_dists_and_trans_by_unit_arma(y),
##   byunit_arma2 = fast_dists_and_trans_by_unit_arma2(y),
##   byunit_par = fast_dists_and_trans_by_unit_arma2_par(y, threads = 7),
##   min_iterations = 100, max_iterations = 1000, check = FALSE, filter_gc = FALSE ## max_iterations to 1000 for below n=1000
## )
## ## Seems like dists_and_trans is much faster with n=5000 (and also with n=100, n=1000, etc..)
## # But it is using a ridiculous amount of memory (382MB versus 354K).
## ## The by unit approaches are so so much better given memory (for n=10000, like 1+GB versus 700 K)
## res1 %>% arrange(median)
# >  res1 (n=100)
##   expression         min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory              time               gc
##   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>              <list>             <list>
## 1 main            29.7µs   46.3µs    21709.  176.27KB        0  1000     0     46.1ms <NULL> <Rprofmem [14 × 3]> <bench_tm [1,000]> <tibble [1,000 × 3]>
## 2 byunit_par     136.4µs  139.5µs     7086.    6.63KB        0  1000     0    141.1ms <NULL> <Rprofmem [6 × 3]>  <bench_tm [1,000]> <tibble [1,000 × 3]>
## 3 byunit_arma2   149.4µs  152.3µs     6505.    6.63KB        0  1000     0    153.7ms <NULL> <Rprofmem [6 × 3]>  <bench_tm [1,000]> <tibble [1,000 × 3]>
## 4 byunit_arma    181.3µs  182.6µs     5416.    6.63KB        0  1000     0    184.6ms <NULL> <Rprofmem [6 × 3]>  <bench_tm [1,000]> <tibble [1,000 × 3]>
## 5 fastdandtrans  227.7µs  270.9µs     3417.    6.63KB        0  1000     0    292.7ms <NULL> <Rprofmem [6 × 3]>  <bench_tm [1,000]> <tibble [1,000 × 3]>
# >  res1 (n=1000)
##   expression         min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory              time             gc
##   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>              <list>           <list>
## 1 main            1.91ms   2.16ms     399.     15.3MB     12.0   200     6   500.93ms <NULL> <Rprofmem [14 × 3]> <bench_tm [200]> <tibble [200 × 3]>
## 2 byunit_par     10.37ms  10.48ms      95.4    41.8KB      0     100     0      1.05s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 3 byunit_arma2   10.51ms  10.63ms      94.0    41.8KB      0     100     0      1.06s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 4 byunit_arma    16.56ms   16.9ms      59.1    41.8KB      0     100     0      1.69s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 5 fastdandtrans  23.97ms  26.89ms      35.9    41.8KB      0     100     0      2.79s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
# (n=5000)
##   expression         min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory              time             gc
##   <bch:expr>    <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>              <list>           <list>
## 1 main              52ms   55.7ms     16.4      382MB     16.6   100   101      6.09s <NULL> <Rprofmem [14 × 3]> <bench_tm [100]> <tibble [100 × 3]>
## 2 byunit_par       271ms  274.3ms      3.63     198KB      0     100     0     27.57s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 3 byunit_arma2     274ms    287ms      3.15     198KB      0     100     0     31.75s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 4 byunit_arma      438ms  439.1ms      2.25     198KB      0     100     0     44.46s <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>
## 5 fastdandtrans    636ms  818.2ms      1.22     198KB      0     100     0      1.36m <NULL> <Rprofmem [6 × 3]>  <bench_tm [100]> <tibble [100 × 3]>

## bench1 <- microbenchmark::microbenchmark(
##   main = dists_and_trans(y),
##   fastdandtrans = fast_dists_and_trans(y),
##   byunit_arma = fast_dists_and_trans_by_unit_arma(y),
##   byunit_arma2 = fast_dists_and_trans_by_unit_arma2(y),
##   byunit_par = fast_dists_and_trans_by_unit_arma2_par(y, threads = 7),
##   times = 100
## )
## bench1
