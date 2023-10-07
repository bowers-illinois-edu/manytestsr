# Test and develop functions to produce p-values.
# Focusing especially on test statistics that will not allow positive and negative effects to cancel out
context("Performance of  P-value Functions")

## The next lines are for use when creating the tests.
## Change interactive<-FALSE for production
interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  library(testthat)
  load_all() ## use  this during debugging
}
setDTthreads(1)

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

  ## Get the names of the transformations of the outcome
  dists <- with(idat[bF == 1], dists_and_trans(Ynull, Z))

  idat[, names(dists) := dists_and_trans(Ynull, Z), by = bF]

  b1[, names(dists) := dists_and_trans(Ynull, Z)]

  Yvers <- c("Ynull", names(dists))
  bigfmla_txt <- paste(paste(Yvers, collapse = "+"), "~ZF|bF", sep = "")
  bigfmla <- as.formula(bigfmla_txt)

  ## maximum of the absolute values of the standardized linear statistic can be more powerful than quadratic
  ## so we use that. Eventually allow people to choose quadratic if they'd like, or even do something fancier.
  it1a <- independence_test(bigfmla, data = idat, teststat = "maximum",alternative="greater")
  pvalue(it1a)
  expect_gt(pvalue(it1a), .05)

  ## Exploring how coin works
  statistic(it1a, type = "linear")
  expectation(it1a)
  corTZ <- zapsmall(cov2cor(covariance(it1a)))
  ## This is the quadratic form
  (statistic(it1a, "linear") - expectation(it1a)) %*% solve(covariance(it1a)) %*% t(statistic(it1a, "linear") - expectation(it1a))
  statistic(it1a)
  ## This is the max form
  max(abs((statistic(it1a, "linear") - expectation(it1a))/diag(covariance(it1a))))

  b1fmla_txt <- paste(paste(Yvers, collapse = "+"), "~ZF", sep = "")
  b1fmla <- as.formula(b1fmla_txt)
  b1ta <- independence_test(b1fmla, data = b1, teststat = "maximum", alternative="greater")
  pvalue(b1ta)
  expect_gt(pvalue(b1ta), .05)

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
    thepsbig[1, i] <- pvalue(independence_test(newbigfmla, data = idat, teststat = "maximum", alternative="greater"))
    for (j in 1:length(Yvers)) {
      newfmla <- as.formula(paste(Yvers[j], "~newZ|bF", sep = ""))
      thepsbig[j + 1, i] <- pvalue(independence_test(newfmla, data = idat, teststat = "maximum", alternative="greater"))
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
    thepsb1[1, i] <- pvalue(independence_test(newb1fmla, data = b1, teststat = "maximum"))
    for (j in 1:length(Yvers)) {
      newfmla <- as.formula(paste(Yvers[j], "~newZ", sep = ""))
      thepsb1[j + 1, i] <- pvalue(independence_test(newfmla, data = b1, teststat = "maximum"))
    }
  }
  fp_b1 <- apply(thepsb1, 1, function(x) {
    mean(x < .05)
  })
  rownames(thepsb1) <- c("Omnibus", Yvers)
  expect_lt(max(fp_b1) - 2 * sqrt((.05 * .95) / 1000), .05)

})

## Observed Effects by Block
idat[, list(
  canceling = mean(Y[Z == 1] - Y[Z == 0]),
  null = mean(Ynull[Z == 1] - Ynull[Z == 0]),
  homog = mean(Yhomog[Z == 1] - Yhomog[Z == 0]),
  normb = mean(Ynormb[Z == 1] - Ynormb[Z == 0])
), by = b]


test_that("pIndepDist works as expected", {
  res_cancel_idat <- pIndepDist(dat = idat, fmla = Y ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_cancel_b1 <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_null_idat <- pIndepDist(dat = idat, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_null_b1 <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_homog_idat <- pIndepDist(dat = idat, fmla = Yhomog ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_homog_b1 <- pIndepDist(dat = b1, fmla = Yhomog ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_normb_idat <- pIndepDist(dat = idat, fmla = Ynormb ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  res_normb_b1 <- pIndepDist(dat = b1, fmla = Ynormb ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore", adaptive_dist_function
  =FALSE)
  expect_lt(res_cancel_idat, .05)
  expect_lt(res_cancel_b1, .05)
  expect_gt(res_null_idat, .05)
  expect_gt(res_null_b1, .05)
  expect_lt(res_homog_idat, .05)
  expect_lt(res_homog_b1, .05)
  expect_lt(res_normb_idat, .05)
  expect_lt(res_normb_b1, .05)
})


 load(file=here::here("tests","dpp_dat.rda"))

 pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF, adaptive_dist_function=FALSE)
 pIndepDist(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
 pWilcox(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
 pWilcox(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
 pOneway(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
 pOneway(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)


independence_test(R01TMCRET~trtF | blockF, data=dpp_dat)
independence_test(rank(R01TMCRET)~trtF | blockF, data=dpp_dat)
ot1 <- oneway_test(R01TMCRET~trtF | blockF, data=dpp_dat)
wt1 <- wilcox_test(R01TMCRET~trtF | blockF, data=dpp_dat)
it1 <- independence_test(R01TMCRET+rank(R01TMCRET)~trtF | blockF, data=dpp_dat)

statistic(ot1)
statistic(wt1)
statistic(it1,"standardized",partial=FALSE)
statistic(it1,"standardized",partial=TRUE)
covariance(it1)
variance(ot1)
variance(wt1)
cov2cor(covariance(it1))
sigma1 <- covariance(it1)
it1_stats <- statistic(it1,"standardized",partial=FALSE)
library(mvtnorm)
dist <- Rfast::rmvnorm(n=1000,mu=rep(0,2),sigma=sigma1)
ps <- matrix(nrow=nrow(dist),ncol=ncol(dist))
ps[,1] <- pnorm(dist[,1],mean=0,sd=sqrt(variance(ot1)))
ps[,2] <- pnorm(dist[,2],mean=0,sd=sqrt(variance(wt1)))
minp <- apply(ps,1,min)
alpha <- 0.05
threshold <- quantile(minp, alpha)
threshold ## a bit higher than pvalue(wt1) but not so high as pvalue(it1)



set.seed(1235)
ot1dist <- rperm(ot1,n=1000)
set.seed(1235)
wt1dist <- rperm(wt1,n=1000)

## from Alex

library(mvtnorm)
library(randomizr)
# Helper functions
do_t_test <- function(Y, Z){
  t.test(Y[Z==1], Y[Z==0])$p.value
}
permute_treat <- function(){
  treatment_sim <- complete_ra(n, m=n/2)
  ps_sim <- apply(outcomes, 2, do_t_test, Z = treatment_sim)
  return(ps_sim)
}
threshold_finder<- function(threshold){
  mean(apply(many_ps, 2, x <- function(x) sum(x <= threshold) > 0 ))
}
# Set a seed
set.seed(343)
# Generate correlated outcomes
# Outcomes are unrelated to treatment
# All null hypotheses are true
n <- 1000
k <- 100; r <- .7; s <- 1
sigma <- matrix(s*r, k,k)
diag(sigma) <- s
outcomes <- rmvnorm(n=n, mean=rep(0, k), sigma=sigma)
# Complete Random Assignment
treatment <- complete_ra(n, m=n/2)
# Conduct k hypothesis tests
p_obs <- apply(outcomes, 2, do_t_test, Z = treatment)
# Simulate under the sharp null
many_ps <- replicate(1000, permute_treat(), simplify = TRUE)
# Obtain the Type I error rate for a series of thresholds
thresholds <- seq(0, 0.05, length.out = 1000)
type_I_rate <- sapply(thresholds, threshold_finder)
# Find the largest threshold that yields an alpha type I error rate
target_p_value <- thresholds[max(which(type_I_rate <=0.05))]
# Apply target p_value to observed p_values
sig_simulated <- p_obs <= target_p_value
# Compare to raw p-values
sig <- p_obs <= 0.05
sig
p_obs[sig]


# L tests
COV <- cov2cor(covariance(it1))
L <- ncol(COV)

# For a range of quantiles (z-values in the context of a standard normal distribution)
z_values <- qnorm(seq(0.001, 0.999, by=0.001))

get_prob_max_z_below <- function(z, COV) {
  lower <- rep(-Inf, L)
  upper <- rep(z, L)
  prob <- pmvnorm(lower = lower, upper = upper, sigma = COV)
  return(prob)
}

prob_max_z_below_values <- sapply(z_values, function(z) get_prob_max_z_below(z, COV))

# Obtain the significance level for p_min
alpha <- 0.05
adjusted_z <- approx(prob_max_z_below_values, z_values, xout = 1 - alpha)$y

# Convert z-value to p-value
adjusted_p_value <- 1 - pnorm(adjusted_z)
print(adjusted_p_value)


dat <- copy(dpp_dat)
outcome_names <- c( "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY")
dat[, c( "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY"):=dists_and_trans(R01TMCRET),by=blockF]

independence_test(R01TMCRET+mndist+mndistRank0+mndistRank0+maddist+maddistRank0+maxdist+zscoreY+rankY~trtF | blockF, data=dat,
		  teststat="quadratic")

idatnew <- copy(idat)
idatnew[, c( "mndist", "mndistRank0", "maddist", "maddistRank0", "maxdist", "maxdistRank0", "zscoreY", "rankY"):=dists_and_trans(Y),by=bF]
pIndepDist(dat = idatnew, fmla = Y~ ZF| bF, adaptive_dist_function=FALSE)

independence_test(Y+mndist+mndistRank0+mndistRank0+maddist+maddistRank0+maxdist+zscoreY+rankY~ZF | bF, data=idatnew,
		  teststat="quadratic")


test_fn <- function(indx,thedat,fmla){
	fmlatxt <- paste(paste(outcomenms[indx],collapse="+"),"~trtF|blockF",collapse="")
	fmla <- as.formula(fmlatxt)
	t1 <- independence_test(fmla,data=thedat,teststat="",distribution=asymptotic())
	pvalue(t1)[[1]]
}

all_but_one <- combn(7,6)
res2 <- apply(all_but_one,2,function(i){ test_fn(indx=i) })
names(res2) <- apply(all_but_one,2,function(idx){ setdiff(1:7,idx) })
res2



test_that("Ordinary tests do not reject when effects cancel across blocks", {
		  expect_gt(pvalue(oneway_test(Y ~ ZF | bF, data = idat)), .05)
		  expect_gt(pvalue(wilcox_test(Y ~ ZF | bF, data = idat)), .05)
		  })

test_that("Ordinary tests do reject when effects are large within blocks even if they cancel across blocks", {
		  expect_lt(max(idat[, list(p = pvalue(wilcox_test(Y ~ ZF, data = .SD))), by = bF]$p), .05)
		  })

test_that("passing a block factor to a p-value function with one block gives the same results as  omitting it", {
		  res_cancel_b1_block <- pIndepDist(dat = b1, fmla = Y ~ ZF | bF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore")
		  res_cancel_b1_noblock <- pIndepDist(dat = b1, fmla = Y ~ ZF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore")
		  expect_equal(res_cancel_b1_block, res_cancel_b1_block)
		  res_null_b1_block <- pIndepDist(dat = b1, fmla = Ynull ~ ZF | bF, distfn = dists_and_trans, , ncpu = 4, parallel = "multicore")
		  res_null_b1_noblock <- pIndepDist(dat = b1, fmla = Ynull ~ ZF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore")
		  expect_equal(res_null_b1_block, res_null_b1_block)
		  })
