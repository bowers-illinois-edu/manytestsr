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
  tmp3 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y, Z = 1)
  tmp5 <- fast_dists_by_unit_arma2_par(y, Z = 1, threads = numcores)

  test_length_fn(tmp1, tmp2)
  test_length_fn(tmp1, tmp3)
  test_length_fn(tmp1, tmp4)
  test_length_fn(tmp1, tmp5)

  test_contents_fn(tmp1, tmp2)
  test_contents_fn(tmp1, tmp3)
  test_contents_fn(tmp1, tmp4)
  test_contents_fn(tmp1, tmp5)

  y <- rchisq(1000, df = 1)
  tmp1 <- dists_and_trans(y)
  tmp2 <- fast_dists_and_trans(y, Z = 1)
  tmp3 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)
  tmp4 <- fast_dists_and_trans_by_unit_arma2(y, Z = 1)
  tmp5 <- fast_dists_by_unit_arma2_par(y, Z = 1, threads = numcores)

  test_length_fn(tmp1, tmp2)
  test_length_fn(tmp1, tmp3)
  test_length_fn(tmp1, tmp4)
  test_length_fn(tmp1, tmp5)

  test_contents_fn(tmp1, tmp2)
  test_contents_fn(tmp1, tmp3)
  test_contents_fn(tmp1, tmp4)
  test_contents_fn(tmp1, tmp5)
})

## ## ### Check timing and core usage
## cores <- parallel::detectCores(logical = FALSE)
## y <- rchisq(100,df=1)
## res1 <- bench::mark(main=dists_and_trans(y),
##     fastdandtrans=fast_dists_and_trans(y, Z = 1),
##     byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##     byunit_arma2=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##     byunit_par=fast_dists_by_unit_arma2_par(y, Z = 1, threads=7),
##     min_iterations=10, max_iterations=100,check=FALSE,filter_gc=FALSE)
## ## Seems like dists_and_trans is much faster with n=5000 (and also with n=100, n=1000, etc..) But it is using a ridiculous amount of memory (382MB versus 354K). The by unit approaches are so so much better given memory (for n=10000, like 1+GB versus 700 K)
## res1
##
## bench1 <- microbenchmark::microbenchmark(main=dists_and_trans(y),
##     fastdandtrans=fast_dists_and_trans(y, Z = 1),
##     byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##     byunit_arma2=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##     byunit_par=fast_dists_by_unit_arma2_par(y, Z = 1, threads=7),
## 	times = 100)
## bench1

## load(file=here::here("Analysis","dpp_dat.rda"))
##
## pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pIndepDist(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
## pWilcox(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pWilcox(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
## pOneway(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
## pOneway(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
##
## ## dpp_dat[,tmpy:=rchisq(.N,df=1)+trt*rnorm(.N),by=blockF]
## ##dpp_dat[,tmpy:=rnorm(.N)+trt*rnorm(.N),by=blockF]
## dpp_dat[,tmpy:=R01TMCRET]
## tmpdat <- copy(dpp_dat)
## tmpdat[,   c("mndist", "mndistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY"):=dists_and_trans(tmpy),by=blockF]
##
## tsmat <- as.matrix(tmpdat[,c("tmpy","mndist", "mndistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")])
## cor(tsmat)
##
## outcomenms <- c("tmpy","mndist", "mndistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")
##
## test_fn <- function(indx){
## fmlatxt <- paste(paste(outcomenms[indx],collapse="+"),"~trtF|blockF",collapse="")
## fmla <- as.formula(fmlatxt)
## t1 <- independence_test(fmla,data=tmpdat,teststat="maximum",distribution=asymptotic())
## pvalue(t1)[[1]]
## }
##
## test_fn(indx=1:length(outcomenms))
## ## Very little penalty for including highly correlated variables.
## test_fn(indx=c(1,1,1:7))
##
## res1 <- sapply(1:length(outcomenms),function(i){ test_fn(indx=i) })
## names(res1) <- outcomenms
## res1
##
## all_but_one <- combn(7,6)
## res2 <- apply(all_but_one,2,function(i){ test_fn(indx=i) })
## names(res2) <- apply(all_but_one,2,function(idx){ setdiff(1:7,idx) })
## res2

## y <- sample(ypop,20000)
## res2 <- bench::mark(#main=dists_and_trans(y),
##     #arma=fast_dists_and_trans(y, Z = 1),
##     byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##     byunit_arma2=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##     byunit_par=fast_dists_by_unit_arma2_par(y, Z = 1, threads=7),
##     min_iterations=2, max_iterations=100,check=FALSE,filter_gc=FALSE)
## res2
