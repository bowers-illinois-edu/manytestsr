## Decide on when to switch among the three different distance matrix functions making the input to independence_test

library(bench)
library(testthat)
devtools::load_all()
setDTthreads(1)

set.seed(123456)
ypop <- rnorm(100000)
## Try to ensure that the values in N=10 are the same 10 (the first 10) as in the N=100, etc.


test_length_fn <- function(obj1,obj2){
    expect_equal(length(obj1),length(obj2))
}

test_contents_fn <- function(obj1,obj2){
    expect_equal(obj1[[1]], obj2[[1]])
    expect_equal(obj1[[2]], obj2[[2]])
    expect_equal(obj1[[3]], obj2[[3]])
    expect_equal(obj1[[4]], obj2[[4]])
    expect_equal(obj1[[5]], obj2[[5]])
    expect_equal(obj1[[6]], obj2[[6]])
    expect_equal(obj1[[7]], obj2[[7]])
    expect_equal(obj1[[8]], obj2[[8]])
}

numcores <- parallel::detectCores(logical=FALSE)
numcores <- 16

test_that("All algorithmns give the same answers",{
y <- 1:5
tmp1 <- dists_and_trans(y)
tmp2 <- fast_dists_and_trans(y, Z = 1)
tmp3 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)
tmp4 <- fast_dists_and_trans_by_unit_arma2(y, Z = 1)
tmp5 <- fast_dists_by_unit_arma2_par(y,Z=1,threads=numcores)
test_length_fn(tmp1,tmp2)
test_length_fn(tmp1,tmp3)
test_length_fn(tmp1,tmp4)
test_length_fn(tmp1,tmp5)
test_contents_fn(tmp1,tmp2)
test_contents_fn(tmp1,tmp3)
test_contents_fn(tmp1,tmp4)
test_contents_fn(tmp1,tmp5)
})

## for more ideas see https://stackoverflow.com/questions/44703599/loop-in-rcpp-slower-than-those-in-r

## y <- sample(ypop,1000)
## n1k_osx_times <- bench::mark(main=dists_and_trans(y),
##         arma=fast_dists_and_trans(y, Z = 1),
##         byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##         byunit_arma=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
##         byunit_arma=fast_dists_by_unit_arma2_par(y, Z = 1,threads=numcores),
##         min_iterations=10, max_iterations=1000,check=FALSE,filter_gc=FALSE)
## n1k_osx_times

dist_timings_small <- press(N=unique(c(seq(10,100,10),seq(100,1000,100))), #,seq(1000,10000,1000),20000),
    {set.seed(12345)
    y <- sample(ypop,N)
    bench::mark(main=dists_and_trans(y),
        arma=fast_dists_and_trans(y, Z = 1),
        byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
        byunit_arma2=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
        byunit_arma2_par=fast_dists_by_unit_arma2_par(y, Z = 1, threads=numcores),
        min_iterations=100, max_iterations=1000,check=FALSE,filter_gc=FALSE)
    })

### UNCOMMENT OUT THESE NEXT TWO FOR OS X LAPTOP TIMINGS
## dist_timings_small_osxlaptop <- dist_timings_small
## save(dist_timings_small_osxlaptop,file="dist_timings_small_osxlaptop.rda")
#### UNCOMMENT OUT THESE NEXT TWO FOR KEELING TIMINGS
## dist_timings_linuxkeeling_small <- dist_timings_small
## save(dist_timings_linuxkeeling,file="dist_timings_linuxkeeling.rda")
dist_timings_small_osxpro <- dist_timings_small
save(dist_timings_small_osxpro,file="dist_timings_small_osxpro.rda")

dist_timings_large <- press(N=unique(c(seq(1000,10000,1000),20000)),
    {set.seed(12345)
    y <- sample(ypop,N)
    bench::mark(main=dists_and_trans(y),
        arma=fast_dists_and_trans(y, Z = 1),
        byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
        byunit_arma2=fast_dists_and_trans_by_unit_arma2(y, Z = 1),
        byunit_arma2_par=fast_dists_by_unit_arma2_par(y, Z = 1, threads=numcores),
        min_iterations=2, max_iterations=1000,check=FALSE,filter_gc=FALSE)
    })

### UNCOMMENT OUT THESE NEXT TWO FOR OS X LAPTOP TIMINGS
## dist_timings_large_osxlaptop <- dist_timings_large
## save(dist_timings_large_osxlaptop,file="dist_timings_large_osxlaptop.rda")
#### UNCOMMENT OUT THESE NEXT TWO FOR KEELING TIMINGS
## dist_timings_large_linuxkeeling<- dist_timings_large
## save(dist_timings_large_linuxkeeling,file="dist_timings_large_linuxkeeling.rda")
dist_timings_large_osxpro <- dist_timings_large
save(dist_timings_large_osxpro,file="dist_timings_large_osxpro.rda")




## Something about fast_dists_and_trans_by_unit is just destroying memory. The byunit_arma is the right approach for big stuff.

##dist_timings_osx2 <- press(N=c(15000,20000),
##    {set.seed(12345)
##    y <- sample(ypop,N)
##    bench::mark(main=dists_and_trans(y),
##        arma=fast_dists_and_trans(y, Z = 1),
##        byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##        ##byunit=fast_dists_and_trans_by_unit(y, Z = 1), ## this would crash the machine
##        min_iterations=10, max_iterations=1000,check=FALSE,filter_gc=FALSE)
##    })
##save(dist_timings_osx2,file="dist_timings_osx2.rda")


## dist_timings_osx3 <- press(N=c(30000,50000,100000),
##     {set.seed(12345)
##     y <- sample(ypop,N)
##     bench::mark(arma=fast_dists_and_trans(y, Z = 1),
##         byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
##         min_iterations=3, max_iterations=1000,check=FALSE,filter_gc=FALSE)
##     })
## save(dist_timings_osx3,file="dist_timings_osx3.rda")

## this next takes too long
## set.seed(12345)
## N <- 800000
## y <- rnorm(N) ##sample(ypop,N)
## n800k_unit <- bench_time( fast_dists_and_trans_by_unit_arma(y, Z = 1) )
## n800k_arma <- bench_time( fast_dists_and_trans(y, Z = 1) )

## pdf(file="dist_timings_osxlaptop.pdf")
## ggplot2::autoplot(dist_timings_osxlaptop)
## dev.off()
