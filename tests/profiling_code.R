## Profiling the code to see if we have any big slowdowns.
## Not really a test. Something to be ignored.
##
##

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- TRUE
if (interactive) {
  library(here)
  library(data.table)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  library(Rfast)
  library(microbenchmark)
  library(bench)
  load_all() ## use  this during debugging
}

load("~/Documents/PROJECTS/manytests-paper/Data/mtwrkdat.rda")

icols <- c("YContNorm", "y0ContNorm", "block", "blockF", "trtF", "trt", "covscluster", "covsplits")
bcols <- c("block", "blockF", "hwt", "truemndiffb", "covscluster", "covsplits")

setkey(mtbdat, block)
setkey(mtwrkdat, block)

set.seed(12345)
blocksB1000 <- sample(mtbdat$block, 1000)
bdatB1000 <- droplevels(mtbdat[.(blocksB1000), ..bcols])
setkey(bdatB1000, block)

idatB1000 <- droplevels(mtwrkdat[.(blocksB1000), ..icols, drop = TRUE])
setkey(idatB1000, block)

## Make relatively sparse effects with average effect of .25sd.
idatB1000$y1 <- create_effects(
  idat = idatB1000, ybase = "y0ContNorm", blockid = "blockF", tau_fn = tau_norm_covariate_outliers,
  tau_size = .25, prop_blocks_0 = .9, covariate = "covscluster"
)
idatB1000[, truetau := (y1 - y0ContNorm)]
truetau_B1000 <- idatB1000[, .(truetau = mean(truetau)), by = blockF]
idatB1000[, Y := trt * y1 + (1 - trt) * y0ContNorm]
idatB1000[, sometau := truetau != 0]
setkey(truetau_B1000, blockF)
setkey(bdatB1000, blockF)
bdatB1000 <- merge(bdatB1000, truetau_B1000)
bdatB1000[, truemndiffb := truetau] ## use the info just created here
## Sometau=TRUE if the true mean difference is not 0
bdatB1000 <- bdatB1000[, sometau := truemndiffb != 0]
table(bdatB1000$sometau)
ate_B1000 <- idatB1000[, .(ate_b = mean(y1 - y0ContNorm)), by = blockF]
setkey(ate_B1000, blockF)
identical(ate_B1000$ate, truetau_B1000$truetau)
stopifnot(unique(round(abs(ate_B1000$ate - truetau_B1000$truetau), 13)) == 0)
setkey(bdatB1000, blockF)

rm(mtbdat)
rm(mtwrkdat)

Rprof(line.profiling = TRUE)
siu_cluster <- findBlocks(
  idat = idatB1000, bdat = bdatB1000, pfn = pIndepDist, splitfn = splitCluster, thealpha = .05,
  fmla = Y ~ trtF | blockF, parallel = "no", trace = TRUE, splitby = "covscluster", stop_splitby_constant = TRUE
)
Rprof(NULL)
summaryRprof(lines = "both")

Rprof(line.profiling = TRUE)
block_tests_fdr <- adjust_block_tests(idat = idatB1000, bdat = bdatB1000, blockid = "block", pfn = pIndepDist, p_adj_method = "fdr", fmla = Y ~ trtF, copydts = TRUE)
Rprof(NULL)
summaryRprof(lines = "both")


set.seed(12345)
xblah <- as.numeric(seq(1.0, 1000.0, length = 2000))
x <- as.numeric(1:10)

### Faster euclidean distances
## First make sure I get it
dist2 <- sqrt(outer(xblah, xblah, FUN = function(x, y) {
  (x - y)^2
}))
testdist <- Dist(xblah)
stopifnot(all.equal(dist2, testdist))

fn0 <- function() { Rfast::vecdist(xblah) }
fn1 <- function() { Dist(xblah) }
fn4 <- function() { vecdist2(xblah) }
fn5 <- function() { vecdist_arma(xblah) }
fn6 <- function() {
  mat <- as.matrix(dist(xblah, diag = TRUE, upper = TRUE))
  dimnames(mat) <- NULL
  return(mat)
}
fn7 <- function() { vecdist3(xblah) }

all.equal(fn6(), fn0())
all.equal(fn6(), fn1())
all.equal(fn6(), fn4())
all.equal(fn6(), fn5())
all.equal(fn6(), fn7())
## microbenchmark(vecdist=fn0(),Dist=fn1(),eigendist=fn2(),fastDist=fn3(),vecdist2=fn4(),vecdist_arma=fn5(),base=fn6(),times=1000)

microbenchmark::microbenchmark(vecdist = fn0(), Dist = fn1(), vecdist2 = fn4(), vecdist3 = fn7(), vecdist_arma = fn5(), base = fn6(), times = 10)

## Another version of vecdist3 but in R
tmpfn <- function(x) {
  C <- -2 * (x %*% t(x))
  An <- matrix(x, ncol = 1)^2
  one <- sweep(C, 2, An, "+")
  two <- sweep(one, 1, An, "+")
  sqrt(two)
}
all.equal(fn1(),tmpfn(xblah))

fn5mat <- vecdist_arma(x)
fn1mat <- Dist(x)

## use manytestsr:::vecdist_arma hmmm or vecdist2 when x is long

## Fast row means
N <- 5000
set.seed(12345)
blah <- matrix(rnorm(N), N / 4, 4)
y <- matrix(blah[, 1], ncol = 1)
y[2:3] <- y[1]

# all.equal(fastrowMeans(blah),fastrowMeans2(blah))
# all.equal(fastrowMeans_Eigen(blah),rowMeans(blah))
all.equal(fastrowMeans(blah), matrix(rowMeans(blah), ncol = 1))
microbenchmark(fastrowMeans(blah), rowMeans(blah), times = 1000)

## use manytestsr::fastrowMeans

### Fast row Mads
tmp <- as.vector(fastrowMads(blah))
tmp0 <- as.vector(fastrowMads2(blah))
tmp2 <- rowMads(blah)
tmp3 <- apply(blah, 1, mad)

all.equal(tmp, tmp2)
all.equal(tmp0, tmp2)
all.equal(tmp2, tmp3)

microbenchmark::microbenchmark(fastrowMads(blah), fastrowMads2(blah), rowMads(blah), apply(blah, 1, mad))

## use manytestsr::fastrowMads2

### Fast row Maxsa
N <- 1000
blah <- matrix(rnorm(N), N / 4, 4)
tmp <- as.vector(fastrowMaxs(blah))
tmp2 <- rowMaxs(blah, value = TRUE)
tmp3 <- apply(blah, 1, max)
tmp4 <- fastrowMaxs2(blah)

all.equal(tmp, tmp2)
all.equal(tmp2, tmp3)
all.equal(tmp2, tmp4)

microbenchmark::microbenchmark(fastrowMaxs(blah), fastrowMaxs2(blah), rowMaxs(blah, value = TRUE), apply(blah, 1, max), times = 1000)

## use manytestsr::fastrowMaxs2


### Fast covariance matrix
all.equal(fastcova(blah), cov(blah))
all.equal(fastcova(blah), cova(blah))
microbenchmark(fastcova(blah), cova(blah), cov(blah))

## Use manytestsr:::fastcova

### Fast variance


tmp1 <- var(y)
tmp2 <- Var(y)
tmp3 <- fastVar(y)
all.equal(tmp1, tmp2)
all.equal(tmp1, tmp3)
microbenchmark(fastVar(y), var(y), Var(y), times = 1000)

## Use manytestsr:::fastVar

### Fast Mahalanobis distance (Which is just a z score) for a vector (not a matrix)
tmp1 <- mahala(y, mu = fastMean(y), sigma = fastcova(y))
tmp2 <- mahalanobis(y, center = fastMean(y), cov = fastcova(y))
tmp4 <- zscore_vec(y)
tmp5 <- as.vector((y - mean(y)) / sd(y))^2
fn6 <- function(y) {
  as.vector(((y - fastMean(y))^2 / Var(y)))
}

all.equal(tmp1, tmp2)
all.equal(tmp1, as.vector(tmp4))
all.equal(tmp1, tmp5)
all.equal(tmp1, fn6(y))
cor(tmp1, tmp5)
cor(tmp1, tmp2)
cor(tmp1, tmp4)

microbenchmark(zscore_vec(y), fn6(y), mahala(y, mu = fastMean(y), sigma = fastcova(y)))

### Ranks

### Make ties
tmp1 <- rank(y, ties.method = "average")
tmp2 <- Rank(y)
tmp4 <- avg_rank_arma(y)
tmp6 <- avg_rank(y)
## tmp5 <- Rank_mean2(y)
all.equal(tmp1, tmp2)
all.equal(tmp1,tmp6)
all.equal(tmp1, as.vector(tmp4))

y <- rnorm(2000)
microbenchmark::microbenchmark(avg_rank_arma(y), avg_rank(y), Rank(y), rank(y, ties.method = "average"), times = 1000)

## avg_rank_arma works. Rank_mean2 is doing something else.

## Faster dists_and_trans
y <- as.numeric(seq(1.0,10.0,length=100))

vecd1 <- Rfast::vecdist(y)
vecd2 <- vecdist_arma(y)
vecd3 <- vecdist2(y)
vecd4 <- vecdist3(y)

all.equal(vecd1, vecd2)
all.equal(vecd1, vecd3)
all.equal(vecd1, vecd4)

tmp1 <- dists_and_trans(y)
tmp2 <- fast_dists_and_trans(y, Z = 1)
tmp3 <- fast_dists_and_trans_by_unit(y, Z = 1)

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

# http://adv-r.had.co.nz/C-interface.html
# https://thecoatlessprofessor.com/programming/cpp/unofficial-rcpp-api-documentation/#iterators

## small y
y <- as.numeric(seq(1.0,100.0))
microbenchmark::microbenchmark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1), times = 100)

## larger y
y <- as.numeric(seq(1.0,1000.0))
microbenchmark::microbenchmark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1), times = 100)

y <- as.numeric(seq(1.0,2000.0))
microbenchmark::microbenchmark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1), times = 100)

## This next the unit by unit approach is better
y <- as.numeric(seq(1.0,10000.0))
microbenchmark::microbenchmark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1), times = 10)

## Can't ruin this next
y <- as.numeric(seq(1.0,50000.0))
microbenchmark::microbenchmark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1), times = 1)

## Try bench::mark 
library(bench)

y <- rnorm(N)
all.equal(dists_and_trans(y), fast_dists_and_trans(y, Z = 1),check.attributes = FALSE)

blah <- press(N=c(100,seq(1000,50000,1000)),
    {y <- rnorm(N)
    bench::mark(dists_and_trans(y), fast_dists_and_trans(y, Z = 1), fast_dists_and_trans_by_unit(y, Z = 1),max_iterations=2,check=FALSE)
    })


## Here we cannot use any of the vecdist functions
ybig <- as.numeric(seq(1.0,10000.0))

## All these fail with 50000 rows
vecd1big <- Rfast::vecdist(ybig)
vecd2big <- vecdist_arma(ybig) ## with latest compiler flag this works but is slow
vecd3big <- vecdist2(ybig)
vecd4big <- vecdist3(ybig)

microbenchmark::microbenchmark( fast_dists_and_trans_by_unit(ybig, Z = 1),fast_dists_and_trans(y, Z = 1), times = 1)

tmpbig0 <- dists_and_trans(ybig) ## this fails

system.time(
tmpbig2 <- fast_dists_and_trans_by_unit(ybig, Z = 1)
)
## 182.999  39.181 222.864  
system.time(
tmpbig1 <- fast_dists_and_trans(ybig, Z = 1)
)

