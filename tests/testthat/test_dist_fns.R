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
setDTthreads(1)

## Make a test vector that is skewed
set.seed(123455)
xblah <- rchisq(100,df=1) + rnorm(100)

test_that("Some distance matrix creation functions produce the same results",{
dist2 <- sqrt(outer(xblah, xblah, FUN = function(x, y) {
  (x - y)^2
}))
testdist <- Dist(xblah)
expect_equal(dist2, testdist)
})

byhand_fn <- function(){ sqrt(outer(xblah, xblah, FUN = function(x, y) { (x - y)^2 })) }
Dist_fn <- function() { Rfast::Dist(xblah,vector=FALSE) }
vecdist_fn <- function() { Rfast::vecdist(xblah) }
#fn1 <- function() { Dist(xblah,square=TRUE) }
manhattan_fn <- function() { Dist(xblah,method="manhattan") }
vecdist2_fn <- function() { vecdist2(xblah) }
vecdist_arma_fn <- function() { vecdist_arma(xblah) }
basedist_fn <- function() {
  mat <- as.matrix(dist(xblah, diag = TRUE, upper = TRUE))
  dimnames(mat) <- NULL
  return(mat)
}
manhattan2_fn <- function(){ manhattan_dist(xblah) }
euc_dist_arma1_fn <- function(){ euc_dist_arma1(xblah) }

test_that("All distance matrix creation functions produce the same results",{
expect_equal(byhand_fn(),Dist_fn())
expect_equal(byhand_fn(),vecdist_fn())
expect_equal(byhand_fn(),manhattan_fn())
expect_equal(byhand_fn(),vecdist2_fn())
expect_equal(byhand_fn(),vecdist_arma_fn())
expect_equal(byhand_fn(),basedist_fn())
expect_equal(byhand_fn(),manhattan2_fn())
expect_equal(byhand_fn(),euc_dist_arma1_fn())
})

bench1 <- bench::mark(byhand=byhand_fn(),vecdist=vecdist_fn(),Dist=Dist_fn(),manhattan=manhattan_fn(),vecdist2=vecdist2_fn(),
	       vecdist_arma=vecdist_arma_fn(),basdist=basedist_fn(),mah2=manhattan2_fn(),euc2=euc_dist_arma1_fn(),
	       min_iterations=100, max_iterations=1000,check=TRUE,filter_gc=TRUE)
bench1 %>% arrange(median)
# # A tibble: 9 × 13
# expression        min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result
# <bch:expr>   <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>
#   1 vecdist        3.65µs     10µs   104315.    86.8KB   209.     998     2     9.57ms <dbl[…]>
#   2 Dist           9.02µs   11.2µs    85724.    88.3KB   172.     998     2    11.64ms <dbl[…]>
#   3 vecdist2      14.43µs   22.7µs    43580.    87.3KB    87.3    998     2     22.9ms <dbl[…]>
#   4 byhand           31µs   43.4µs    21528.     235KB   108.     995     5    46.22ms <dbl[…]>
#   5 manhattan     48.38µs   53.9µs    18568.    88.1KB    18.6    999     1     53.8ms <dbl[…]>
#   6 vecdist_arma  62.03µs   68.4µs    14804.    87.3KB    29.7    998     2    67.41ms <dbl[…]>
#   7 basdist       69.45µs   81.7µs    11970.   444.5KB   121.     990    10     82.7ms <dbl[…]>
#   8 mah2         136.28µs  149.6µs     6618.    87.3KB     6.62   999     1   150.94ms <dbl[…]>
#   9 euc2         160.43µs  165.1µs     5881.    87.3KB    11.8    998     2    169.7ms <dbl[…]>

## bench2 <- microbenchmark::microbenchmark(byhand=byhand_fn(),vecdist=vecdist_fn(),Dist=Dist_fn(),manhattan=manhattan_fn(),vecdist2=vecdist2_fn(),
## 	       vecdist_arma=vecdist_arma_fn(),basdist=basedist_fn(),mah2=manhattan2_fn(),euc2=euc_dist_arma1_fn(),
## 	       times=1000,check="equal")
##
## bench2df <- summary(bench2) %>% arrange(median)
# n=100, vecdist,Dist, vecdist2, byhand
# n=1000, same
# n=10000, same

best4a <- bench1 %>% arrange(median) %>%  filter(median <= median[4])
# # A tibble: 4 × 13
# expression     min  median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result   memory
# <bch:expr> <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list>   <list>
#   1 vecdist     3.53µs  6.11µs   159496.    86.8KB    320.    998     2     6.26ms <dbl[…]> <Rprofmem>
#   2 Dist        9.14µs 11.77µs    84778.    88.3KB     84.9   999     1    11.78ms <dbl[…]> <Rprofmem>
#   3 vecdist2   16.24µs 24.23µs    41881.    87.3KB     41.9   999     1    23.85ms <dbl[…]> <Rprofmem>
#   4 byhand     30.05µs 39.28µs    23689.     235KB     71.3   997     3    42.09ms <dbl[…]> <Rprofmem>
#   # ℹ 2 more variables: time <list>, gc <list>
# best4b <- bench2df %>% arrange(median) %>% filter(median <= median[4])
#      expr    min     lq      mean median     uq      max neval     cld
# 1  vecdist  3.567  5.043  6.397435  6.560  7.134   46.289  1000   c
# 2     Dist  8.979 10.865 34.917937 11.685 12.341 4947.716  1000  bc
# 3 vecdist2 20.623 22.878 28.810167 24.026 25.010 4600.241  1000   cd
# 4   byhand 30.135 37.433 62.677848 39.565 41.615 4863.830  1000 ab


test_that("Best distance algorithms are the four that are best on development machine",{
		  expect_equal(as.character(best4a$expression),c("vecdist","Dist","vecdist2","byhand"))
	       })
