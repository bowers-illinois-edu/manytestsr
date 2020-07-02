# Test and develop functions to produce p-values.
# Focusing especially on test statistics that will not allow positive and negative effects to cancel out
context("Performance of  P-value Functions")

## The next lines are for use when creating the tests. Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  ## library(here)
  ## library(data.table)
  ## library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  load_all() ## use  this during debugging
}

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

test_that("Basics of the omnibus test works as expected in regards False Positive Rate for single tests", {
  ### Ensuring no effects
  ## idat[, Z := sample(Z), by = bF]
  ## idat[, ZF := factor(Z, labels = c(0, 1))]

  dists <- with(idat[bF == 1], dists_and_trans(Ynull, Z))

  idat[, c(names(dists), "rankY") := c(
    dists_and_trans(Ynull, Z),
    list(Rfast::Rank(Ynull)) # ,
    # list(Rfast::mahala(matrix(Ynull, ncol = 1), mu = mean(Ynull), sigma = cov(matrix(Ynull, ncol = 1))))
  ), by = bF]
  ## idat[,rankmhY:=Rfast::mahala(matrix(rankY,ncol=1), mu=mean(rankY),sigma=cov(matrix(rankY,ncol=1))),by=bF]

  b1[, c(names(dists), "rankY") := c(
    dists_and_trans(Ynull, Z),
    list(Rfast::Rank(Ynull)) # ,
    # list(Rfast::mahala(matrix(Ynull, ncol = 1), mu = mean(Ynull), sigma = cov(matrix(Ynull, ncol = 1))))
  )]

  Yvers <- c("Ynull", names(dists), "rankY")
  bigfmla_txt <- paste(paste(Yvers, collapse = "+"), "~Z|bF", sep = "")
  bigfmla <- as.formula(bigfmla_txt)

  it1a <- independence_test(bigfmla, data = idat, teststat = "quadratic")
  pvalue(it1a)
  expect_gt(pvalue(it1a), .05)

  b1fmla_txt <- paste(paste(Yvers, collapse = "+"), "~Z", sep = "")
  b1fmla <- as.formula(b1fmla_txt)
  b1ta <- independence_test(b1fmla, data = b1, teststat = "quadratic")
  pvalue(b1ta)
  expect_gt(pvalue(b1ta), .05)

  ## statistic(it1a, type = "linear")
  ## expectation(it1a)
  ## corTZ <- zapsmall(cov2cor(covariance(it1a)))
  ## corYnull <- zapsmall(cor(idat[, .SD, .SDcols = Yvers]))
  ## all.equal(dimnames(corTZ), dimnames(corYnull))
  ## corTZ - corYnull
  ### Done above
  ## Exclude rankmhY  since cor with mndistRank0==1.
  ## bigfmla <- update(bigfmla,.-rankmhY~.)
  ## Yvers <- c("Ynull",names(dists),"rankY","mhY")
  ## bigfmla_txt <- paste(paste(Yvers,collapse="+"),"~Z|bF",sep="")
  ## bigfmla <- as.formula(bigfmla_txt)
  ## it1b <- independence_test(bigfmla,data=idat,teststat="quadratic")
  ## pvalue(it1b)

  ## Assess false positive rate of the individual test (not FWER) on the big dataset
  ## This takes a while right now
  set.seed(12345)
  nsims <- 1000
  thepsbig <- matrix(NA, ncol = nsims, nrow = length(Yvers) + 1)
  newbigfmla_txt <- paste(paste(Yvers, collapse = "+"), "~newZ|bF", sep = "")
  newbigfmla <- as.formula(newbigfmla_txt)
  for (i in 1:nsims) {
    idat[, newZ := sample(Z), by = bF]
    thepsbig[1, i] <- pvalue(independence_test(newbigfmla, data = idat, teststat = "quadratic"))
    for (j in 1:length(Yvers)) {
      newfmla <- as.formula(paste(Yvers[j], "~newZ|bF", sep = ""))
      thepsbig[j + 1, i] <- pvalue(independence_test(newfmla, data = idat, teststat = "quadratic"))
    }
  }
  fp_idat <- apply(thepsbig, 1, function(x) {
    mean(x < .05)
  })
  rownames(thepsbig) <- c("Omnibus", Yvers)
  expect_lt(max(fp_idat) - 2 * sqrt((.05 * .95) / 1000), .05)

  thepsb1 <- matrix(NA, ncol = nsims, nrow = length(Yvers) + 1)
  newb1fmla_txt <- paste(paste(Yvers, collapse = "+"), "~newZ", sep = "")
  newb1fmla <- as.formula(newb1fmla_txt)
  for (i in 1:nsims) {
    b1[, newZ := sample(Z), by = bF]
    thepsb1[1, i] <- pvalue(independence_test(newb1fmla, data = b1, teststat = "quadratic"))
    for (j in 1:length(Yvers)) {
      newfmla <- as.formula(paste(Yvers[j], "~newZ", sep = ""))
      thepsb1[j + 1, i] <- pvalue(independence_test(newfmla, data = b1, teststat = "quadratic"))
    }
  }
  fp_b1 <- apply(thepsb1, 1, function(x) {
    mean(x < .05)
  })
  rownames(thepsb1) <- c("Omnibus", Yvers)
  expect_lt(max(fp_b1) - 2 * sqrt((.05 * .95) / 1000), .05)
})

test_that("pIndepDist works as expected", {
  res_cancel_idat <- pIndepDist(dat = idat, fmla = Y ~ ZF | bF, distfn = dists_and_trans)
  res_cancel_b1 <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = dists_and_trans)
  res_null_idat <- pIndepDist(dat = idat, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans)
  res_null_b1 <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans)
  res_homog_idat <- pIndepDist(dat = idat, fmla = Yhomog ~ ZF | bF, distfn = dists_and_trans)
  res_homog_b1 <- pIndepDist(dat = b1, fmla = Yhomog ~ ZF | bF, distfn = dists_and_trans)
  res_normb_idat <- pIndepDist(dat = idat, fmla = Ynormb ~ ZF | bF, distfn = dists_and_trans)
  res_normb_b1 <- pIndepDist(dat = b1, fmla = Ynormb ~ ZF | bF, distfn = dists_and_trans)
  expect_lt(res_cancel_idat, .05)
  expect_lt(res_cancel_b1, .05)
  expect_gt(res_null_idat, .05)
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
  expect_lt(max(idat[, list(p = pvalue(wilcox_test(Y ~ ZF, data = .SD))), by = bF]$p), .05)
})

test_that("passing a block factor to a p-value function with one block gives the same results as  omitting it", {
  res_cancel_b1_block <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = dists_and_trans)
  res_cancel_b1_noblock <- pIndepDist(dat = b1, fmla = Y ~ ZF, distfn = dists_and_trans)
  expect_equal(res_cancel_b1_block, res_cancel_b1_block)
  res_null_b1_block <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans)
  res_null_b1_noblock <- pIndepDist(dat = b1, fmla = Ynull ~ ZF, distfn = dists_and_trans)
  expect_equal(res_null_b1_block, res_null_b1_block)
})
