## Make some testing data for use in the tests in this directory

library(randomizr)
library(data.table)
setDTthreads(1)
set.seed(12345)

## For now, create a very simple canceling out arrangement and simple design setup
idat <- data.table(i = 1:1000, b = rep(c(1:10), length = 1000))
bdat <- data.table(b = 1:10)
setkey(bdat, b)
setkey(idat, b)
## Make vb1, vb2, vb3 all nest within each other like states, counties, districts (with the individual blocks nested within all three)
## 2 levels at highest level
bdat[, vb1 := rep(c(1, 2), length.out = .N)]
## 2 levels within each state
bdat[, vb2 := rep(1:2, length.out = .N), by = vb1]
bdat[, vbnest := interaction(vb1, vb2, lex.order = TRUE, drop = TRUE)]
ftable(bdat$vb1, bdat$vb2)
## Now merge the block data onto the individual level data
idat <- bdat[idat]
set.seed(12345)
idat[, vi1 := round(runif(.N, min = 0, max = 10)), by = b]
idat[, y0 := vi1 + vb1 + rnorm(100), by = b]
idat[, tau := ifelse(b <= 5, -5 * sd(y0), 5 * sd(y0))]
idat[, tauhomog := sd(y0) * 5]
idat[, taunormb := rnorm(.N, mean = (sd(y0) + sd(y0) / b), sd = sd(y0) / 2), by = b]
idat[, y1 := y0 + tau]
idat[, y1homog := y0 + tauhomog]
idat[, y1normb := y0 + taunormb]
stopifnot(all.equal(mean(idat$y1 - idat$y0), 0))
idat[, Z := complete_ra(N = 100), by = b]
testra <- idat[, sum(Z), by = b]
stopifnot(all(testra$V1 == 50))
idat[, Y := Z * y1 + (1 - Z) * y0]
## Now ensure that the no effects case is really no effects
idat[, y0null := resid(lm(y0 ~ Z, data = .SD)), by = b]
idat[, y1null := y0null]
idat[, Ynull := Z * y1null + (1 - Z) * y0null]
stopifnot(all.equal(idat$Ynull, idat$y0null))
idat[, Yhomog := Z * y1homog + (1 - Z) * y0]
idat[, Ynormb := Z * y1normb + (1 - Z) * y0]
idat$ZF <- factor(idat$Z)
idat$bF <- factor(idat$b)
idat$gF <- interaction(idat$bF, idat$ZF)
## Created aligned versions of the  variables
idat[, Ymd := Y - mean(Y), by = b]
idat[, Zmd := Z - mean(Z), by = b]

## Make some blocks have nothing going on so the splitting doesn't go all the way to the end
set.seed(12345)
idat[bF %in% c(7, 8), Ynormb := Ynull + rnorm(.N)]

## Observed Effects by Block
idat[, list(
  canceling = mean(Y[Z == 1] - Y[Z == 0]),
  null = mean(Ynull[Z == 1] - Ynull[Z == 0]),
  homog = mean(Yhomog[Z == 1] - Yhomog[Z == 0]),
  normb = mean(Ynormb[Z == 1] - Ynormb[Z == 0])
), by = b]

##  Make a smaller dataset
b1n6 <- droplevels(idat[b %in% c(1, 6), ])
setorder(b1n6, b)
b1n6
b1 <- droplevels(b1n6[b == 1, ])

bdat_tmp <- idat[, .(nb = .N, nt = sum(Z), nc = sum(1 - Z), pb = mean(Z)), by = b]
bdat_tmp[, hwt := (2 * (nc * nt) / (nc + nt))]
setkey(bdat_tmp, b)
bdat <- bdat_tmp[bdat]
bdat$bF <- factor(bdat$b)

setkey(idat, bF)
setkey(bdat, bF)

#### Now make data with unequal sized blocks, unequal taus, and more blocks
bdat2 <- bdat[, .(nb = round(seq(1, 20, length = 20), 2), bF2 = factor(1:20))]
setkey(bdat2, "bF2")
idat[, b2 := {
  sel <- rep(c(0, 1), .N / 2)
  sel * b + (1 - sel) * (10 + b)
}, by = bF]
idat[, bF2 := factor(b2)]
setkey(idat, "bF2")
## Make a dataset with unequal sized blocks
idat2 <- bdat2[idat]
idat3 <- idat2[, sample(.N, size = (.N * nb) / 10, replace = TRUE), by = bF2]
## Now giving it the same blocking names as idat for ease with testing
idat3$bF <- idat3$bF2
setkey(idat3, bF)
idat3 <- idat3[, b := as.numeric(as.character(bF2))]
idat3[, c("v1", "v2", "v3", "v4") := lapply(1:4, function(i) {
  rnorm(.N, mean = .N / i, sd = 1 / i)
}), by = bF2]
idat3[, y0 := v1 + rnorm(.N), by = bF2]
idat3[, tau := ifelse(b <= 10, -5 * sd(y0), 5 * sd(y0))]
idat3[, tauhomog := 5 * sd(y0)]
## Taunorm_inc increases the effect with the block number and also size of the block.
## We should be increasingly likely to detect effects for blocks with higher numbers.
idat3[, taunorm_inc := rnorm(.N, mean = b / 10 * sd(y0), sd = .5), by = bF2]
idat3[, taunorm_inc := ifelse(taunorm_inc < 0, 0, taunorm_inc)]
## taunorm_dec has effect size decreasing in block size
idat3[, taunorm_dec := rnorm(.N, mean = 4 / b * sd(y0), sd = .5), by = bF2]
idat3[, taunorm_dec := ifelse(taunorm_dec < 0, 0, taunorm_dec)]
idat3[, sd(y0), by = b]
## tauv2 has no systematic relationship between block number or block size and effect size
idat3[, tauv2 := ifelse(v2 < quantile(v2, sample(c(.9, .1), size = 1)), sd(y0) * (taunorm_inc) + sd(y0) * runif(.N, min = .25, max = 2), 0), by = bF2]
idat3[, lapply(.SD, mean), .SDcols = grep("tau|nb|y0", names(idat3), value = TRUE), by = b]
idat3[, y1 := y0 + tau]
idat3[, y1null := y0]
idat3[, y1homog := y0 + tauhomog]
set.seed(1235) ## set y1norm_inc=y0 for 1/4 of the blocks
nullbF <- sample(unique(levels(idat3$bF)), size = length(unique(idat3$bF)) / 4)
idat3[, y1norm_inc := ifelse(bF %in% nullbF, y0, y0 + taunorm_inc)]
idat3[, y1norm_dec := ifelse(bF %in% nullbF, y0, y0 + taunorm_dec)]
idat3[, y1tauv2 := y0 + tauv2]
idat3[, nb := .N, by = bF2]
set.seed(1235)
tmp <- idat3[, .(nb = .N), by = bF2]
mb <- round(tmp$nb * runif(20, min = .2, max = .8))
idat3[, Z := block_ra(blocks = bF2, block_m = mb)]
idat3[, Y := Z * y1 + (1 - Z) * y0]
idat3[, Ynull := Z * y1null + (1 - Z) * y0]
idat3[, Yhomog := Z * y1homog + (1 - Z) * y0]
idat3[, Ynorm_inc := Z * y1norm_inc + (1 - Z) * y0]
idat3[, Ynorm_dec := Z * y1norm_dec + (1 - Z) * y0]
idat3[, Ytauv2 := Z * y1tauv2 + (1 - Z) * y0]
idat3$ZF <- factor(idat3$Z)
idat3[, .(.N, mean(y1tauv2 - y0)), by = bF]
bdat3 <- idat3[, .(
  nb = .N, nt = sum(Z), nc = sum(1 - Z), pb = mean(Z),
  v1 = mean(v1), v2 = mean(v2), v3 = mean(v3), v4 = mean(v4),
  ate_tau = mean(y1 - y0),
  ate_null = mean(y1null - y0),
  ate_homog = mean(y1homog - y0),
  ate_norm_inc = mean(y1norm_inc - y0),
  ate_norm_dec = mean(y1norm_dec - y0),
  ate_tauv2 = mean(y1tauv2 - y0)
), by = bF]
bdat3[, hwt := (2 * (nc * nt) / (nc + nt))]

setkey(idat, bF)
setkey(bdat, bF)
setkey(idat3, bF)
setkey(bdat3, bF)
