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
  load_all() ## use  this during debugging
}

load("~/Documents/PROJECTS/manytests-paper/data/mtwrkdat.rda")

icols <- c("YContNorm", "y0ContNorm", "block", "blockF", "trtF", "trt", "covscluster", "covsplits")
bcols <- c("block", "blockF", "hwt", "truemndiffb", "covscluster", "covsplits")

setkey(mtbdat, block)
setkey(mtwrkdat, block)

set.seed(12345)
blocksB1000 <- sample(mtbdat$block,1000)
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


Rprof(line.profiling=TRUE)
siu_cluster <- findBlocks(
  idat = idatB1000, bdat = bdatB1000, pfn = pIndepDist, splitfn = splitCluster, thealpha = .05,
  fmla = Y ~ trtF | blockF, parallel = "no", trace=TRUE, splitby="covscluster", stop_splitby_constant = TRUE
)
Rprof(NULL)

Rprof(line.profiling=TRUE)
block_tests_fdr <- adjust_block_tests(idat=idatB1000,bdat=bdatB1000,blockid="block",pfn=pIndepDist,p_adj_method="fdr",fmla=Y~trtF,copydts=TRUE)
Rprof(NULL)

rdist <- function(x){
  G <- t(x) %*% x
  D <- diag(G) + t(diag(G)) - 2*G
  return(D)
}

m1 <- c(0.68,-0.211,0.566,0.597,0.823)
m2 <- c(-0.605,-0.33,0.536,-0.444,0.108)
Dist(cbind(m1,m2))

x <- 1:10
G <- x %*% t(x)
xsq <- diag(G)*2
Dfin <- 2*G

dist2 <- sqrt(outer(x,x,FUN=function(x,y){ (x-y)^2 }))
testdist <- Dist(x)
all.equal(dist2,testdist)

# set.seed(12345)
xblah <- rnorm(1000)
fn0 <- function(){ vecdist(xblah) }
fn1 <- function(){ Dist(xblah) }
all.equal(fn1(),fn0())
fn2 <- function(){ eigenDist(xblah) }
all.equal(fn1(),fn2())
fn3 <- function(){ calcPWD1(matrix(xblah,ncol=1)) }
all.equal(fn1(),fn3())
microbenchmark(fn0(),fn1(),fn2(),fn3(),times=100)
#
# Rprof(line.profiling=TRUE)
# blah <- dists_and_trans(x=xblah)
# Rprof(NULL)
# summaryRprof(lines="both")

N <- 200
blah  <- matrix(rnorm(N),N/4,4)
all.equal(fastrowMeans(blah),row_means(blah))
all.equal(fastrowMeans(blah),matrix(rowMeans(blah),ncol=1))
microbenchmark(fastrowMeans(blah),row_means(blah),rowMeans(blah))


all.equal(fastcova(blah),cov(blah))
all.equal(fastcova(blah),cova(blah))
microbenchmark(fastcova(blah),cova(blah),cov(blah))

tmp <- as.vector(fastrowMads(blah))
tmp0 <- as.vector(fastrowMads2(blah))
tmp2 <- rowMads(blah)
tmp3 <- apply(blah,1,mad)

all.equal(tmp,tmp2)
all.equal(tmp0,tmp2)
all.equal(tmp2,tmp3)

microbenchmark(fastrowMads(blah),fastrowMads2(blah),rowMads(blah),apply(blah,1,mad))

