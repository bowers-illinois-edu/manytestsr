# Test and develop functions to produce p-values.
# Focusing especially on test statistics that will not allow positive and negative effects to cancel out
context("Performance of  P-value Functions")

library(data.table)
library(testthat)
library(here)
library(coin)
source(here("tests/testthat", "make_test_data.R"))
devtools::load_all()


set.seed(12345)
dat <- data.frame(y0 = sample(1:7), z = sample(c(1, 1, 1, 0, 0, 0, 0)))
dat$y1 <- dat$y0 + 1
dat$Y <- with(dat, z * y1 + (1 - z) * y1)
dat <- dat[order(dat$z), ]
dat$z
table(dat$z)

dy <- dist(dat$Y)
dymat <- as.matrix(dy)
dymat[lower.tri(dymat)] <- 0
dymat2 <- as.matrix(dy)

library(energy)
ed0a <- edist(dat$Y, sizes = c(4, 3))
ed0b <- eqdist.e(dat$Y, sizes = c(4, 3))
ed0c <- ksample.e(dat$Y, sizes = c(4, 3))
ed1 <- Rfast::edist(dat$Y[dat$z == 1], dat$Y[dat$z == 0])
all.equal(ed0a[[1]], ed0b[[1]])
all.equal(ed0a[[1]], ed0c[[1]])
ed1

#### Now my version
nT <- sum(dat$z)
nC <- sum(1 - dat$z)

Tix <- which(dat$z == 1)
Cix <- which(dat$z == 0)

dymat2T <- dymat2[Tix, Tix]
dymat2C <- dymat2[Cix, Cix]
dymat2TC <- dymat2[Tix, Cix]
totT <- sum(dymat2T)
totC <- sum(dymat2C)
totTC <- sum(dymat2TC)
avgT <- totT / (nT * nT) ## m22
avgC <- totC / (nC * nC) ## m11
avgTC <- totTC / (nT * nC) ## m12
thew <- (nT * nC) / (nT + nC) ## w
e_mine <- thew * ((avgTC + avgTC) - (avgC + avgT))
all.equal(ed0a[[1]], e_mine)

## For speed we might use the functions from Rfast
totTb <- Rfast::total.dista(dat$Y[dat$z == 1], dat$Y[dat$z == 1])
totCb <- Rfast::total.dista(dat$Y[dat$z == 0], dat$Y[dat$z == 0])
totTCb <- Rfast::total.dista(dat$Y[dat$z == 1], dat$Y[dat$z == 0])
all.equal(totCb, totC)
all.equal(totTb, totT)
all.equal(totTCb, totTC)

### Now, can we get the same number or a number proportional to it by first working at the individual level?

TCi <- colSums(dymat2TC)
all.equal(sum(TCi), totTC)

TCiC <- colSums(dymat2TC)
TCiT <- rowSums(dymat2TC)
all.equal(sum(TCiC), sum(TCiT))
all.equal(sum(TCiC), totTC)

Ti <- colSums(dymat2T) ## the T and C mats will be square so col or row is same
Ci <- colSums(dymat2C) ## the T and C mats will be square so col or row is same
all.equal(sum(Ti), totT)
all.equal(sum(Ci), totC)

wi <- ((nT * nC) / (nT + nC)) / nrow(dat)

e_mine2 <- thew * ((sum(TCiC) / (nT * nC) - sum(Ci) / (nC * nC)) + (sum(TCiT) / (nT * nC) - sum(Ti) / (nT * nT)))
all.equal(ed0a[[1]], e_mine2)

z <- dat$z
dat$id <- 1:nrow(dat)
names(z) <- dat$id

eCi <- ((TCiC / (nT * nC)) - (Ci / (nC * nC))) ## assuming same order
eTi <- ((TCiT / (nT * nC)) - (Ti / (nT * nT))) ## assuming same order

ei <- c(eCi, eTi)

e_mine3 <- thew * sum(ei)
all.equal(e_mine3, ed0a[[1]])

e_mine4 <- sum(ei * thew)

## So e per person is something like this

ei2 <- c(eCi * thew, eTi * thew)
all.equal(sum(ei2), ed0a[[1]])

newei <- edisti(dat$Y, dat$z)


## also try just mh distance so that we can calculate it without reference to Z

## mvnfast::maha
## Rfast::mahala
## stats::mahalanobis

ytmp <- rnorm(1000)
Ymat <- matrix(ytmp, ncol = 1)
centers <- mean(ytmp)
blah <- microbenchmark(stats::mahalanobis(Ymat, center = centers, cov = cov(Ymat)),
  Rfast::mahala(Ymat, mu = centers, sigma = cov(Ymat)),
  mvnfast::maha(Ymat, mu = centers, sigma = cov(Ymat)),
  times = 1000
)


## Or something like mean or median of pairwise distances

# Also rank based mahalanobis

##  We have Y where the effects  are gigantic within  block and canceling
## Ynull where the sharp  null is  true
## Yhomog where all blocks have additive  constant effect of  5
## And Ynormb where blocks have hetergeneous  effects following a normal  dist within block
## We want a procedure that, for Y, rejects the null both within any given block as well as across blocks.
##  The simple tests fail to achieve this goal
# oneway_test(Y~ZF|bF,data=idat) ## what we do not  want
# oneway_test(Y~ZF|bF,data=b1) ## what we want
# wilcox_test(Y~ZF|bF,data=idat) ## what we do not  want
# wilcox_test(Y~ZF|bF,data=b1) ## what we want

# b1 is block level
# idat is individual level


### ## Neither of these should detect effects
set.seed(12345)
idat[, Z := sample(Z), by = bF]
idat[, ZF := factor(Z, labels = c(0, 1))]

dists <- with(idat[bF == 1], ctrldist(Ynull, Z))

## This next does not work
idat[, c("e_i1", "e_i2", "e_i3", "e_i4", "e_i5", "e_i6", "rankY") := c(
  ctrldist(Ynull, Z),
  list(Rfast::Rank(Ynull))
), by = bF]
it1a <- independence_test(e_i1 + e_i2 + e_i3 + e_i4 + e_i5 + e_i6 + Ynull + rankY ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1a)
statistic(it1a, type = "linear")
expectation(it1a)
covariance(it1a)
cov2cor(covariance(it1a))

it1b <- independence_test(Ynull ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1b)
it1c <- independence_test(e_i1 ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1c)
it1d <- independence_test(e_i2 ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1d)
it1e <- independence_test(e_i3 ~ Z | bF, data = idat, teststat = "quadratic") ## less good
pvalue(it1e)
it1f <- independence_test(e_i4 ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1f)
it1g <- independence_test(e_i5 ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1g)
it1h <- independence_test(e_i6 ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(it1h)

idat[, c("e_Y_i1", "e_Y_i2", "e_Y_i3", "e_Y_i4", "e_Y_i5", "e_Y_i6", "rankY_Y") := c(
  ctrldist(Y, Z),
  list(Rfast::Rank(Y))
), by = bF]
itYa <- independence_test(e_Y_i1 + e_Y_i2 + e_Y_i3 + e_Y_i4 + e_Y_i5 + e_Y_i6 + Y + rankY_Y ~ Z | bF, data = idat, teststat = "quadratic")
pvalue(itYa) # this should be small. the other ones should be big.
statistic(itYa, type = "linear")
expectation(itYa)
covariance(itYa)
cov2cor(covariance(itYa))

b1 <- droplevels(idat[bF == 1, ])
it0a <- independence_test(Ynull + rankY ~ Z, data = b1, teststat = "quadratic")
statistic(it0a)
statistic(it0a, type = "linear")
statistic(it0a, type = "standardized")
b1[, .(sum(Ynull * Z), sum(rankY * Z))]
it0b <- independence_test(Ynull + rankY ~ ZF, data = b1, teststat = "quadratic")
statistic(it0b, type = "linear")
b1[, .(sum(Ynull * (ZF == "0")), sum(rankY * (ZF == "0")))]
statistic(it0b, type = "standardized")
statistic(it0b)
it0c <- independence_test(e_i1 ~ Z, data = b1, teststat = "quadratic", distribution = approximate(nresample = 1000))
support(it0c)
statistic(it0c, type = "linear")


library(Ball)
blah <- bd.test(Ynull ~ Z, data = b1, method = "limit", num.permutations = 1000)

xb0a <- RItools::xBalance(Z ~ Ynull + rankY, data = b1, report = "all")
xb0a
xb0b <- RItools::xBalance(Z ~ Ynull + rankY + e_i1 + e_i2 + e_i3, data = b1, report = "all")
xb0b

setorder(b1, Z)
b0 <- as.data.frame(b1)
b0 <- b0[c(1:5, 95:100), ]
ctrldist(b0$Ynull, b0$Z)
dx <- Dist(b0$Ynull)
dz <- outer(b0$Z, b0$Z, function(x, y) {
  abs(x - y)
})
M <- cbind(b0$Z, (1 - b0$Z))
M_sums <- dx %*% M
Mxy_sums <- (M_sums * (1 - M)) %*% c(1, 1)
Mxx_sums <- (M_sums * M) %*% c(1, 1)

Mij_i <- rowSums(dx)
Mii_i <- rowSums(dx)
blah <- edist(dist(b1$Ynull), sizes = c(50, 50))
blah <- ksample.e(dist(b1$Ynull), sizes = c(50, 50), distance = TRUE)
blah2 <- ksample.e(dist(b1$Ynull), sizes = c(50, 50), distance = TRUE, method = "discoF")
blah3 <- ksample.e(dist(b1$Ynull), sizes = c(50, 50), distance = TRUE, method = "discoB") ## original*2 because uses full sample denom

# xb1 <- RItools::xBalance(Z~rankY+rank_e_i+e_i+Ynull,strata=list(bf=~bF),data=idat,report="all")
# xb1$overall
# expect_equal(pvalue(it1a)[[1]],xb1$overall[,3])

res_null_idat <- pIndepDist(dat = idat, fmla = Ynull ~ Z | bF, distfn = ctrldist)
res_null_idat

it4a <- independence_test(rankY ~ ZF | bF, data = idat, teststat = "quadratic")
it4b <- independence_test(rank_e_i ~ ZF | bF, data = idat, teststat = "quadratic")
it4c <- independence_test(e_i ~ ZF | bF, data = idat, teststat = "quadratic")
it4d <- independence_test(Ynull ~ ZF | bF, data = idat, teststat = "quadratic")
pvalue(it4a)
pvalue(it4b)
pvalue(it4c)
pvalue(it4d)

nsims <- 1000
thepsbig <- matrix(NA, ncol = nsims, nrow = 8)
for (i in 1:nsims) {
  idat[, newZ := sample(Z), by = bF]
  idat[, c("e_i1", "e_i2", "e_i3", "e_i4", "e_i5", "e_i6", "rankY") := c(ctrldist(Ynull, newZ), list(Rfast::Rank(Ynull)), list(Rfast::mahala(Ynull, mu = mean(Ynull), cov = cov(Ynull)))), by = bF]
  idat[, newZF := factor(newZ)]
  thepsbig[1, i] <- pvalue(independence_test(e_i1 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[2, i] <- pvalue(independence_test(e_i2 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[3, i] <- pvalue(independence_test(e_i3 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[4, i] <- pvalue(independence_test(e_i4 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[5, i] <- pvalue(independence_test(e_i5 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[6, i] <- pvalue(independence_test(e_i6 ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[7, i] <- pvalue(independence_test(Ynull ~ newZF | bF, data = idat, teststat = "quadratic"))
  thepsbig[8, i] <- pvalue(independence_test(rankY ~ newZF | bF, data = idat, teststat = "quadratic"))
}
apply(thepsbig, 1, function(x) {
  mean(x < .05)
})

#############
set.seed(12345)
b1$Y2 <- rnorm(nrow(b1), mean = mean(b1$Ynull))

b1[, c(
  "mndist", "mndistRank0", "maddist", "maddistRank0",
  "maxdist", "maxdistRank0", "edist", "rankY", "mh"
) := c(
  ctrldist(Y2, newZ),
  list(edisti(Y2, newZ)),
  list(Rfast::Rank(Y2)),
  list(Rfast::mahala(matrix(Y2, ncol = 1),
    mu = mean(Y2), sigma = cov(matrix(Y2, ncol = 1))
  ))
)]
b1[, rankmh := Rfast::mahala(matrix(rankY, ncol = 1), mu = mean(rankY), sigma = cov(matrix(rankY, ncol = 1)))]

Yvers <- c(
  "mndist", "mndistRank0", "maddist", "maddistRank0",
  "maxdist", "maxdistRank0", "edist", "rankY", "mh", "rankmh", "Y2"
)
zapsmall(cor(b1[, .SD, .SDcols = Yvers]))

nsims <- 1000
theps <- matrix(NA, ncol = nsims, nrow = length(Yvers))
row.names(theps) <- Yvers
for (i in 1:nsims) {
  b1[, newZ := sample(Z)]
  b1[, c(
    "mndist", "mndistRank0", "maddist", "maddistRank0",
    "maxdist", "maxdistRank0", "edist", "rankY", "mh"
  ) := c(
    ctrldist(Y2, newZ),
    list(edisti(Y2, newZ)),
    list(Rfast::Rank(Y2)),
    list(Rfast::mahala(matrix(Y2, ncol = 1),
      mu = mean(Y2), sigma = cov(matrix(Y2, ncol = 1))
    ))
  )]
  b1[, rankmh := Rfast::mahala(matrix(rankY, ncol = 1), mu = mean(rankY), sigma = cov(matrix(rankY, ncol = 1)))]
  b1[, newZF := factor(newZ)]
  for (j in 1:length(Yvers)) {
    newfmla <- reformulate("newZF", response = Yvers[j])
    ## print(newfmla)
    theps[j, i] <- pvalue(independence_test(newfmla, data = b1, teststat = "quadratic"))
  }
}

apply(theps, 1, function(x) {
  mean(x < .05)
})
## Notes: e_i2 and e_i3 produce same p values, also e_i4 and e_i5 produce same
## pvalues

blah <- independence_test(e_i1 ~ newZF, data = b1, teststat = "quadratic", distribution = approximate(nresample = 1000))
statistic(blah)

## This does not because Z is not breaking its relationship directly with the Distances
nsims <- 1000
ps2 <- rep(NA, nsims)
for (i in 1:nsims) {
  ps2[i] <- pIndepDist(fmla = Ynull ~ newZ | bF, dat = idat[, newZ := sample(Z), by = bF], distfn = ctrldist)
}
mean(ps2 < .05)
hist(ps2)


idat$y1test_null <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 0, prop_blocks_0 = 1, covariate = "vb1")
expect_equal(idat$y1test_null, idat$y0)
idat[, Ynull2 := Z * y1test_null + (1 - Z) * y0]
expect_equal(idat$Ynull2, idat$y0)
all.equal(idat$Ynull, idat$Ynull2)

idat$y0_simp <- rnorm(nrow(idat), mean = 0, sd = 1)
idat$y1_simp <- idat$y0_simp
idat$Y_simp <- idat$y0_simp




idat3[, rank_e_i := list(ctrldist(Ynull, Z)), by = bF]
idat3[, rankY := Rfast::Rank(Ynull), by = bF]
it99 <- independence_test(rankY + rank_e_i + Ynull ~ as.factor(Z) | bF, data = idat3, teststat = "quadratic")
xb99 <- RItools::xBalance(Z ~ rankY + rank_e_i + Ynull, strata = list(bf = ~bF), data = idat3, report = "all")
pvalue(it99)
xb99$overall

cors <- idat[,
  {
    thecor <- cor(.SD)
    max(thecor[lower.tri(thecor)])
  },
  by = bF,
  .SDcols = c("rankY", "e_i", "rank_e_i", "Ynull")
]


res_null_idat3 <- pIndepDist(dat = idat3, fmla = Ynull ~ ZF | bF, distfn = ctrldist)


# it2 <- independence_test(rankY+rank_e_i+Y_simp~as.factor(Z),data=droplevels(idat[bF==5,]),teststat="quadratic")
# xb2 <- RItools::xBalance(Z~rankY+rank_e_i+Y_simp,data=droplevels(idat[bF==5,]),report="all")
# pvalue(it2)
# xb2$overall

ps <- idat[, pvalue(independence_test(rankY + rank_e_i + Y_simp ~ as.factor(Z), data = droplevels(.SD), teststat = "quadratic")), by = bF]
psXB <- idat[, RItools::xBalance(Z ~ rankY + rank_e_i + Y_simp, data = droplevels(.SD), report = "all")$overall, by = bF]
cbind(psXB, ps)

cors <- idat[,
  {
    thecor <- cor(.SD)
    max(thecor[lower.tri(thecor)])
  },
  by = bF,
  .SDcols = c("rankY", "rank_e_i", "Y_simp")
]

thecor <- cor(data.frame(idat[bF == 6, .(rankY, rank_e_i, Ynull)]))
max(thecor[lower.tri(thecor, diag = FALSE)])

summary(idat[bF == 6, .(Z, rankY, rank_e_i, Ynull)])
summary(idat[bF == 7, .(Z, rankY, rank_e_i, Ynull)])

ctrldist(idat[bF == 6, Ynull], idat[bF == 6, Z])


test_that("pIndepDist works as expected", {
  res_cancel_idat <- pIndepDist(dat = idat, fmla = Y ~ ZF | bF, distfn = ctrldist)
  res_cancel_b1 <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = ctrldist)
  res_null_idat <- pIndepDist(dat = idat, fmla = Ynull ~ ZF | bF, distfn = ctrldist)
  res_null_idat3 <- pIndepDist(dat = idat3, fmla = Ynull ~ ZF | bF, distfn = ctrldist)
  ## res_null2_idat <- pIndepDist(dat = idat, fmla = Ynull2 ~ ZF | bF, distfn = ctrldist)
  ## res_null_idat_T <- pOneway(dat = idat, fmla = Ynull ~ ZF | bF)
  res_null_b1 <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = ctrldist)
  res_homog_idat <- pIndepDist(dat = idat, fmla = Yhomog ~ ZF | bF, distfn = ctrldist)
  res_homog_b1 <- pIndepDist(dat = b1, fmla = Yhomog ~ ZF | bF, distfn = ctrldist)
  res_normb_idat <- pIndepDist(dat = idat, fmla = Ynormb ~ ZF | bF, distfn = ctrldist)
  res_normb_b1 <- pIndepDist(dat = b1, fmla = Ynormb ~ ZF | bF, distfn = ctrldist)
  expect_lt(res_cancel_idat, .05)
  expect_lt(res_cancel_b1, .05)
  expect_gt(res_null_idat, .05)
  expect_gt(res_null_idat3, .05)
  expect_gt(res_null_b1, .05)
  expect_lt(res_homog_idat, .05)
  expect_lt(res_homog_b1, .05)
  expect_lt(res_normb_idat, .05)
  expect_lt(res_normb_b1, .05)
})

test_that("Ordinary tests do not reject when effects cancel across blocks", {
  expect_gt(pvalue(oneway_test(Y ~ ZF | bF, data = idat)), .05)
  expect_gt(pvalue(wilcox_test(Y ~ ZF | bF, data = idat)), .05)
})

test_that("Ordinary tests do reject when effects are large within blocks even if they cancel across blocks", {
  expect_lt(max(idat[, list(p = pvalue(wilcox_test(Y ~ ZF, data = .SD))), by = b]$p), .05)
})

test_that("passing a block factor to a p-value function with one block gives the same results as  omitting it", {
  res_cancel_b1_block <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = ctrldist)
  res_cancel_b1_noblock <- pIndepDist(dat = b1, fmla = Y ~ ZF, distfn = ctrldist)
  expect_equal(res_cancel_b1_block, res_cancel_b1_block)
  res_null_b1_block <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = ctrldist)
  res_null_b1_noblock <- pIndepDist(dat = b1, fmla = Ynull ~ ZF, distfn = ctrldist)
  expect_equal(res_null_b1_block, res_null_b1_block)
})
