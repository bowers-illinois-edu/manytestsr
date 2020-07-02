## Decide on when to switch among the three different distance matrix functions making the input to independence_test

library(bench)
devtools::load_all()
setDTthreads(1)

set.seed(123456)
ypop <- rnorm(100000)
## Try to ensure that the values in N=10 are the same 10 (the first 10) as in the N=100, etc.
### For os x, crashes on N=20000
dist_timings <- press(N=c(10,100,500,1000,5000,10000,15000),
    {set.seed(12345)
    y <- sample(ypop,N)
    bench::mark(main=dists_and_trans(y),
        arma=fast_dists_and_trans(y, Z = 1),
        byunit=fast_dists_and_trans_by_unit(y, Z = 1),min_iterations=10, max_iterations=1000,check=FALSE,filter_gc=FALSE)
    })
dist_timings_osxlaptop <- dist_timings
save(dist_timings_osxlaptop,file="dist_timings_osxlaptop.rda")


## Where does it crash?
set.seed(12345)
N <- 15000
y <- sample(ypop,N)
n15k_main <- bench_time( dists_and_trans(y) )
n15k_arma <- bench_time( fast_dists_and_trans(y, Z = 1) )
n15k_unit <- bench_time( fast_dists_and_trans_by_unit(y, Z = 1) )

set.seed(12345)
N <- 20000
y <- sample(ypop,N)
n20k_main <- bench_time( dists_and_trans(y) )
n20k_arma <- bench_time( fast_dists_and_trans(y, Z = 1) )
n20k_unit <- bench_time( fast_dists_and_trans_by_unit(y, Z = 1) )

## We also did this on keeling:
load("dist_timings_linuxkeeling.rda")



