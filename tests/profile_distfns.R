## Decide on when to switch among the three different distance matrix functions making the input to independence_test

library(bench)
devtools::load_all()
setDTthreads(1)

set.seed(123456)
ypop <- rnorm(100000)
## Try to ensure that the values in N=10 are the same 10 (the first 10) as in the N=100, etc.


y <- 1:5

tmp1 <- dists_and_trans(y)
tmp2 <- fast_dists_and_trans(y, Z = 1)
tmp3 <- fast_dists_and_trans_by_unit(y, Z = 1)
tmp4 <- fast_dists_and_trans_by_unit_arma(y, Z = 1)

all.equal(length(tmp1), length(tmp2))
all.equal(length(tmp1), length(tmp3))
all.equal(tmp1[[1]], tmp2[[1]])
all.equal(tmp1[[2]], tmp2[[2]])
all.equal(tmp1[[3]], tmp2[[3]])
all.equal(tmp1[[4]], tmp2[[4]])
all.equal(tmp1[[5]], tmp2[[5]])
all.equal(tmp1[[6]], tmp2[[6]])
all.equal(tmp1[[7]], tmp2[[7]])
all.equal(tmp1[[8]], tmp2[[8]])

all.equal(tmp1[[1]], as.vector(tmp3[[1]]))
all.equal(tmp1[[2]], as.vector(tmp3[[2]]))
all.equal(tmp1[[3]], as.vector(tmp3[[3]]))
all.equal(tmp1[[4]], as.vector(tmp3[[4]]))
all.equal(tmp1[[5]], as.vector(tmp3[[5]]))
all.equal(tmp1[[6]], as.vector(tmp3[[6]]))
all.equal(tmp1[[7]], as.vector(tmp3[[7]]))
all.equal(tmp1[[8]], as.vector(tmp3[[8]]))

all.equal(tmp1[[1]], as.vector(tmp4[[1]]))
all.equal(tmp1[[2]], as.vector(tmp4[[2]]))
all.equal(tmp1[[3]], as.vector(tmp4[[3]]))
all.equal(tmp1[[4]], as.vector(tmp4[[4]]))
all.equal(tmp1[[5]], as.vector(tmp4[[5]]))
all.equal(tmp1[[6]], as.vector(tmp4[[6]]))
all.equal(tmp1[[7]], as.vector(tmp4[[7]]))
all.equal(tmp1[[8]], as.vector(tmp4[[8]]))


y <- sample(ypop,1000)
n1k_osx_times <- bench::mark(main=dists_and_trans(y),
        arma=fast_dists_and_trans(y, Z = 1),
        byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
        byunit=fast_dists_and_trans_by_unit(y, Z = 1),
        min_iterations=10, max_iterations=1000,check=FALSE,filter_gc=FALSE)
n1k_osx_times

dist_timings <- press(N=c(10,100,500,1000,5000,10000,15000),
    {set.seed(12345)
    y <- sample(ypop,N)
    bench::mark(main=dists_and_trans(y),
        arma=fast_dists_and_trans(y, Z = 1),
        byunit_arma=fast_dists_and_trans_by_unit_arma(y, Z = 1),
        byunit=fast_dists_and_trans_by_unit(y, Z = 1),
        min_iterations=10, max_iterations=1000,check=FALSE,filter_gc=FALSE)
    })
##dist_timings_osxlaptop <- dist_timings
##save(dist_timings_osxlaptop,file="dist_timings_osxlaptop.rda")
dist_timings_linuxkeeling <- dist_timings
save(dist_timings_linuxkeeling,file="dist_timings_linuxkeeling.rda")

## pdf(file="dist_timings_osxlaptop.pdf")
## ggplot2::autoplot(dist_timings_osxlaptop)
## dev.off()


## We also did this on keeling:
## load("dist_timings_linuxkeeling.rda")

## summary(dist_timings_linuxkeeling)

pdf(file="dist_timings_linuxkeeling.pdf")
ggplot2::autoplot(dist_timings_linuxkeeling)
dev.off()

