# Test and develop functions to produce p-values.
# Focusing especially on test statistics that will not allow positive and negative effects to cancel out
# And on outcomes that people would treat as continuous but which have lots of ties and zeros and skew
context("Performance of  P-value Functions")

interactive <- FALSE
if (interactive) {
  library(here)
  library(data.table)
  library(devtools)
  source(here("tests/testthat", "make_test_data.R"))
  library(testthat)
  library(dplyr)
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  library(parallel)
  library(mvtnorm)
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
  it1a <- independence_test(bigfmla, data = idat, teststat = "quadratic")
  pvalue(it1a)
  expect_gt(pvalue(it1a), .05)
  ## Exploring how coin works
  statistic(it1a, type = "linear")
  coin::expectation(it1a)
  corTZ <- zapsmall(cov2cor(covariance(it1a)))
  ## This is the quadratic form
  (statistic(it1a, "linear") - coin::expectation(it1a)) %*% solve(covariance(it1a)) %*% t(statistic(it1a, "linear") - coin::expectation(it1a))
  statistic(it1a,type="test")
  ## This is the max form
  it1a_max <- independence_test(bigfmla,data=idat,teststat="maximum")
  statistic(it1a,type="test")
  ## Fix this next
  max(abs((statistic(it1a_max, "linear") - coin::expectation(it1a_max)))/sqrt(variance(it1a_max)))
  b1fmla_txt <- paste(paste(Yvers, collapse = "+"), "~ZF", sep = "")
  b1fmla <- as.formula(b1fmla_txt)
  b1ta <- independence_test(b1fmla, data = b1, teststat = "quadratic")
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


## Below is the search for a better set of test statistics:
## For now using:
## outcome_names <- c(theresponse,"mndist","mndistRank0","maxdist","rankY","tanhx")
## See pval_fns.R, dists.R, and code for fast_dists... in the fastfns.cpp file

## Now using more discrete type data typical of policy applications
## Here we have a lot of zeros which creates ties and/or makes measures of spread degenerate --- no variance within certain blocks for example.
## Find the actual test stats
 load(file=here::here("tests","dpp_dat.rda"))


 test_that("Test using data with lots of zeros and ties", {
   dpp_R01 <- pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF, distfn = dists_and_trans, ncpu = 4, parallel = "multicore")
   expect_lt(dpp_R01,.05)
 })

table(dpp_dat$R01TMCRET, exclude = c())
table(dpp_dat$R02TMCRET, exclude = c())
#
# ## ## Ranks more powerful. R02 easier to detect effects.
#  pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF, adaptive_dist_function=FALSE)
#  pIndepDist(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF, adaptive_dist_function=TRUE)
#  pIndepDist(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
#  pWilcox(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
#  pWilcox(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
#  pOneway(dat = dpp_dat, fmla = R01TMCRET ~ trtF | blockF)
#  pOneway(dat = dpp_dat, fmla = R02TMCRET ~ trtF | blockF)
#
# ## Checking the code above: pOneway and pWilcox should be like these
# ot0 <- independence_test(R01TMCRET~trtF | blockF, data=dpp_dat)
# wt0 <- independence_test(rank(R01TMCRET)~trtF | blockF, data=dpp_dat)
#
# ot1 <- oneway_test(R01TMCRET~trtF | blockF, data=dpp_dat)
# wt1 <- wilcox_test(R01TMCRET~trtF | blockF, data=dpp_dat)
#
# stopifnot(all.equal(pvalue(ot0),pvalue(ot1)))
# stopifnot(all.equal(pvalue(wt0),pvalue(wt1)))
#
# ## See where the loss of power is happening
# dat <- copy(dpp_dat)
# dat[,Y:=R01TMCRET]
#
# dists <- with(idat[bF == 1], dists_and_trans(Ynull, Z))
# outcome_names <- c("Y",names(dists))
# rm(dists)
# dat[, (outcome_names[-1]):=dists_and_trans(R01TMCRET),by=blockF]
#
# idatnew <- copy(idat)
# idatnew[, (outcome_names[-1]):=dists_and_trans(Y),by=bF]
#
# ## Compare the combined with one by one.
# ## A function to allow us to Try subsets of the transformed variables on both the more continuous and less continuous data
# ## Is mad and max either irrelevant or a problem?
# ## Can we get away with fewer functions of the outcome and its distances?
# test_fn <- function(indx,thedat,the_outcome_names=outcome_names,raw_outcome_nm,trtnm,blocknm){
#   ## Assume that outcome_names[1] is "Y"
# 	outcome_names0 <- c(raw_outcome_nm,the_outcome_names[-1])
# 	if(length(indx)>1){
# 	fmlatxt <- paste(paste(outcome_names0[indx],collapse="+"),"~",trtnm,"|",blocknm,collapse="")
# 	} else {
# 	fmlatxt <- paste(outcome_names0[indx],"~",trtnm,"|",blocknm,collapse="")
# 	}
# 	fmla <- as.formula(fmlatxt)
# 	t1 <- independence_test(fmla,data=thedat,teststat="quadratic",distribution=asymptotic())
# 	res0 <- as.data.frame(matrix(ncol=length(outcome_names0)+1,nrow=1))
# 	names(res0) <- c(outcome_names0,"p")
# 	res0[,outcome_names0[indx]] <- outcome_names0[indx]
# 	res0$p=pvalue(t1)[[1]]
# 	return(res0)
# }
#
#
# res_dppR01_singles<- sapply(1:length(outcome_names),function(i){ test_fn(indx=i,thedat=dat,raw_outcome_nm="Y",trtnm="trtF",blocknm="blockF")})
# colnames(res_dppR01_singles) <- outcome_names
# single_ps <- unlist(res_dppR01_singles["p",] )
# single_ps
# p.adjust(single_ps,method="holm")
#
# y_fmlatxt <- as.formula(paste(paste(outcome_names,collapse="+"),"~trtF|blockF",collapse=""))
# it_dat_quad <- independence_test(y_fmlatxt,data=dat,teststat="quadratic")
# it_dat_max <- independence_test(y_fmlatxt,data=dat,teststat="maximum")
#
# pvalue(it_dat_quad)
# pvalue(it_dat_max)
#
# ## The basic test statistics and their implied covariances are all the same
# it_dat_quad_stat <- drop(statistic(it_dat_quad,"standardized",partial=TRUE)) ## from array [1,1:9,1:44] to just [1:9,1:44]
# it_dat_max_stat <- drop(statistic(it_dat_max,"standardized",partial=TRUE))
# all.equal(it_dat_max_stat,it_dat_quad_stat)
#
# statistic(it_dat_quad,type="test")
# statistic(it_dat_max,type="test")
#
# thevar <- variance(it_dat_quad)
# thecov <- covariance(it_dat_quad)
# invcov <- zapsmall(covariance(it_dat_quad,invert=TRUE))
# T <- statistic(it_dat_quad,type="linear")
# mu <- coin::expectation(it_dat_quad)
# thecorr <- cov2cor(thecov)
#
#
# ## > thecorr
# ##                     0:Y 0:mndist 0:mndistRank0 0:maddist 0:maddistRank0 0:maxdist 0:maxdistRank0 0:mhdist  0:rankx 0:mnsqrtdist  0:hubmn 0:tanhx
# ## 0:Y             1.00000  0.80618       0.69137  -0.05553       -0.06874  -0.50720       -0.36379   0.6128  0.79585      0.87396  0.80434  0.8772
# ## 0:mndist        0.80618  1.00000       0.75016   0.19360        0.14980  -0.15188       -0.00376   0.9211  0.54384      0.92755  0.98826  0.5893
# ## 0:mndistRank0   0.69137  0.75016       1.00000   0.31795        0.36262  -0.06669        0.10842   0.6932  0.75607      0.67494  0.70962  0.4653
# ## 0:maddist      -0.05553  0.19360       0.31795   1.00000        0.82319   0.54095        0.60591   0.3441 -0.09813     -0.03605  0.14529 -0.3527
# ## 0:maddistRank0 -0.06874  0.14980       0.36262   0.82319        1.00000   0.46541        0.74331   0.2865 -0.13508     -0.05656  0.11170 -0.3181
# ## 0:maxdist      -0.50720 -0.15188      -0.06669   0.54095        0.46541   1.00000        0.77554   0.1899 -0.49103     -0.47476 -0.20998 -0.7750
# ## 0:maxdistRank0 -0.36379 -0.00376       0.10842   0.60591        0.74331   0.77554        1.00000   0.2470 -0.51786     -0.28744 -0.03871 -0.6096
# ## 0:mhdist        0.61275  0.92112       0.69323   0.34415        0.28648   0.18988        0.24702   1.0000  0.36781      0.72871  0.89417  0.2957
# ## 0:rankx         0.79585  0.54384       0.75607  -0.09813       -0.13508  -0.49103       -0.51786   0.3678  1.00000      0.63954  0.52223  0.7269
# ## 0:mnsqrtdist    0.87396  0.92755       0.67494  -0.03605       -0.05656  -0.47476       -0.28744   0.7287  0.63954      1.00000  0.93630  0.8163
# ## 0:hubmn         0.80434  0.98826       0.70962   0.14529        0.11170  -0.20998       -0.03871   0.8942  0.52223      0.93630  1.00000  0.6110
# ## 0:tanhx         0.87718  0.58929       0.46531  -0.35269       -0.31811  -0.77498       -0.60957   0.2957  0.72689      0.81632  0.61104  1.0000
# ##
# myquadt <- as.vector(T - mu) %*% invcov %*% as.vector(T - mu)
# ## This is the max test stat.
# ## The covariances don't matter just the variances.
# statistic(it_dat_max,type="test")
# mymaxt <- max(abs( (T - mu)/ sqrt(thevar) ) )
#
# set.seed(1234)
# pvalue(it_dat_max)
# set.seed(1234)
# pvalue(it_dat_max,method="global",distribution="joint")
# pvalue(it_dat_max,method="global",distribution="marginal",type="Sidak")
# pvalue(it_dat_max,method="single-step",distribution="joint",type="Sidak")
# pvalue(it_dat_max,method="step-down",distribution="joint",type="Sidak")
# pvalue(it_dat_max,method="single-step",distribution="marginal",type="Sidak")
# pvalue(it_dat_max,method="step-down",distribution="marginal",type="Sidak")
# ## Two ways to get the unadjusted p-values
# pvalue(it_dat_max,method="unadjusted")
# single_ps
#
# ## Notice that the multivariate normal evaluation uses some randomness in the default algorithm
# set.seed(12345)
# 1 - coin:::pmvn(lower=-abs(mymaxt),upper=abs(mymaxt),mean=rep(0,ncol(thecorr)), corr=thecorr,conf.int=FALSE)
# set.seed(12345)
# 1- pmvnorm(lower=-abs(mymaxt),upper=abs(mymaxt),mean=rep(0,ncol(thecorr)), corr=thecorr,keepAttr = FALSE)
# ## This next is much slower
# # 1- pmvnorm(lower=-abs(mymaxt),upper=abs(mymaxt),mean=rep(0,ncol(thecorr)), corr=thecorr,keepAttr = FALSE, algorithm="Miwa")
# set.seed(12345)
# pvalue(it_dat_max,distribution="joint",method="global")[1]
#
# library(MASS)
#
# ## Now try without mad or max
# newcorr <- thecorr[c(1:3,8:9),c(1:3,8:9)]
# 1- pmvnorm(lower=-abs(mymaxt),upper=abs(mymaxt),mean=rep(0,ncol(newcorr)), corr=newcorr,keepAttr = FALSE)
#
# drop(statistic(it_dat_quad,"standardized",partial=FALSE))
# drop(statistic(it_dat_max,"standardized",partial=FALSE))
#
# ## Does it matter that some are negative?
# ## If we removed the mads (with low correlations and thus perhaps higher adjustment) would we do a better job?
# cov_it_dat_quad_stat <- covariance(it_dat_quad)
# cov_it_dat_max_stat <- covariance(it_dat_max)
# all.equal(cov_it_dat_max_stat,cov_it_dat_quad_stat)
#
# options(digits=4)
# zapsmall(cov2cor(cov_it_dat_quad_stat))
# zapsmall(cov2cor(cov_it_dat_max_stat))
#
# ## Some blocks have constant outcomes and so all test statistics are NaN
# dat %>% filter(blockF=="Block0092")
# ## Others are small: hard to calc some of this stuff in pairs like mean distances (difference in mean distance between trt and ctrl is 0)
# dat %>% filter(blockF=="Block0109")
# ## Other have little variance
# dat %>% filter(blockF=="Block0110")
#
# ## Try without the mads or maybe maxes. Check below for all possible combinations
#
# num_cols <- 1:(length(outcome_names))
# poss_cols0 <- lapply(num_cols,function(i){ combn(length(outcome_names),num_cols[i],simplify=FALSE) })
# poss_cols <- unlist(poss_cols0,recursive=FALSE)
#
# ncores <- detectCores()-1
# options(mc.cores=ncores)
#
# res_cancel_lst <- mclapply(poss_cols,function(i){ test_fn(indx=i,thedat=idatnew,raw_outcome_nm="Y",trtnm="ZF",blocknm="bF")})
# res_homog_lst<- mclapply(poss_cols,function(i){ test_fn(indx=i,thedat=idatnew,raw_outcome_nm="Yhomog",trtnm="ZF",blocknm="bF")})
# res_normb_lst<- mclapply(poss_cols,function(i){ test_fn(indx=i,thedat=idatnew,raw_outcome_nm="Ynormb",trtnm="ZF",blocknm="bF")})
#
# res_dppR01_lst <- mclapply(poss_cols,function(i){ test_fn(indx=i,thedat=dat,raw_outcome_nm="Y",trtnm="trtF",blocknm="blockF")})
# res_dppR02_lst <- mclapply(poss_cols,function(i){ test_fn(indx=i,thedat=dat,raw_outcome_nm="R02TMCRET",trtnm="trtF",blocknm="blockF")})
#
# res_cancel <- dplyr::bind_rows(res_cancel_lst)
# res_homog <- dplyr::bind_rows(res_homog_lst)
# res_normb <- dplyr::bind_rows(res_normb_lst)
# res_dppR01 <- dplyr::bind_rows(res_dppR01_lst)
# res_dppR02 <- dplyr::bind_rows(res_dppR02_lst)
#
# num_not_na <- function(x){ sum(!is.na(x)) }
#
# ## Check for degenerate statistic? (attr(mat,"r") >= num obs  - num strata from RItools)
#
# ## Lots of ways to reject the null in this case even without mndist
# res_cancel <- res_cancel %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# summary(res_cancel$p)
# table(res_cancel$num_stats)
# ## For res_cancel including everything works great.
# quantile(res_cancel$p,seq(0,1,.1))
# res_cancel_high_power <- res_cancel %>% filter(p <= .05) %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# nrow(res_cancel_high_power)
# res_cancel_high_power %>% summarize(across(one_of(outcome_names),num_not_na))
# res_cancel_high_power %>% head()
#
# res_cancel_low_power <- res_cancel %>% filter(p >= .1) %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# res_cancel_low_power
# res_cancel_low_power %>% summarize(across(one_of(outcome_names),num_not_na),n())
#
# homog_outcome_names <- c("Yhomog",outcome_names[-1])
# res_homog <- res_homog %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(homog_outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# quantile(res_homog$p,seq(0,1,.1))
# res_homog_high_power <- res_homog %>% filter(p <= .05)
# res_homog_low_power <- res_homog %>% filter(p > .05)
# table(res_homog_low_power$num_stats,exclude=c())
# table(res_homog_high_power$num_stats,exclude=c())
#
# res_homog_low_power %>% filter(num_stats>=6)
# res_homog_high_power %>% filter(num_stats>=8)
#
# normb_outcome_names <- c("Ynormb",outcome_names[-1])
# res_normb <- res_normb %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(normb_outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# quantile(res_normb$p,seq(0,1,.1))
# res_normb_high_power <- res_normb %>% filter(p <= .05)
# res_normb_low_power <- res_normb %>% filter(p > .05)
# table(res_normb_low_power$num_stats,exclude=c())
# res_normb_low_power %>% filter(num_stats >= 6)
# res_normb_high_power %>% filter(num_stats>=8)
#
#
# res_dppR01 <- res_dppR01 %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
# ## Few test stats work with this outcome with so many zeros
# quantile(res_dppR01$p,seq(0,1,.1))
# res_dppR01_high_power <- res_dppR01 %>% filter(p <= .05)
# res_dppR01_low_power <- res_dppR01 %>% filter(p > .05)
# ## Only a few test stats seem to work
# table(res_dppR01_high_power$num_stats,exclude=c())
# table(res_dppR01_low_power$num_stats,exclude=c())
# res_dppR01_low_power
# res_dppR01_high_power %>% filter(num_stats>=3)
# res_dppR01_high_power %>% arrange(desc(num_stats),p)
# quantile(res_dppR01_high_power$p,seq(0,1,.1))
# ## The p-values are barely less than .05 with the different sets of three  stats (.04..)
# ## The lower p-values arise with 2 and 1 test stat.
# ## Only one uses the raw outcome:
# res_dppR01_high_power %>% filter(!is.na(Y))
# res_dppR01_high_power %>% filter(!is.na(rankx))
# res_dppR01_high_power %>% filter(!is.na(mndist))
#
# dpp02_outcome_names <- c("R02TMCRET",outcome_names[-1])
# res_dppR02 <- res_dppR02 %>% rowwise() %>%
# 	mutate(num_stats=num_not_na(c_across(one_of(dpp02_outcome_names)))) %>% ungroup() %>% arrange(desc(num_stats),p)
#
# library(tidyverse)
#
# res_dppR01$Y <- ifelse(!is.na(res_dppR01$Y),"R01TMCRET",NA)
# res_dppR02 <- res_dppR02 %>% rename(Y=R02TMCRET)
# res_homog <- res_homog %>% rename(Y=Yhomog)
# res_normb <- res_normb %>% rename(Y=Ynormb)
#
# res <- bind_rows(list(dppR01=res_dppR01,
# 		      dppR02=res_dppR02,
# 		      cancel=res_cancel,
# 		      homog=res_homog,
# 		      normb=res_normb),.id="Outcome")
#
# ## For the simulated datasets, more test statistics equals more power
# ## For the integer and zero inflated datasets, more test statistics is less power --- the rank based tests are clearly superior
# ## For now, choose a collection that works for all --- not as powerful as the best for each outcome, but still powerful.
#
# ## Possible: do 10 univariate tests and choose the minimum p-value after
# ## adjusting for the the correlation (which might be more powerful than using
# ## the linear combination).
#
# g1 <- res %>%
# 	ggplot(aes(x=num_stats,y=p))+
# 	geom_point()+
# 	facet_grid(~Outcome)
# g1
#
#
# g2 <- filter(res, p<.05) %>%
# 	ggplot(aes(x=num_stats,y=p))+
# 	geom_point()+
# 	facet_grid(~Outcome)
# g2
#
# ## Require at least 4 test stats
# ## Many of the p=0 for cancel, so, choose the version with the most stats for them
# res %>% group_by(Outcome) %>% filter(p <= .05) %>% filter(num_stats>=4) %>%
#   filter(p==min(p)) %>% filter(num_stats==max(num_stats))
#
# res %>% group_by(Outcome) %>% filter(p <= .05) %>% filter(num_stats>=5 & num_stats <=8) %>%
# 	filter(p==min(p)) %>% filter(num_stats==max(num_stats))
#
#
# ##    Outcome Y         mndist mndistRank0 maddist maddistRank0 maxdist maxdistRank0 mhdist rankx mnsqrtdist hubmn tanhx       p num_stats
# ##    <chr>   <chr>     <chr>  <chr>       <chr>   <chr>        <chr>   <chr>        <chr>  <chr> <chr>      <chr> <chr>   <dbl>     <int>
# ##  1 dppR01  R01TMCRET NA     NA          NA      NA           maxdist NA           mhdist rankx NA         NA    tanhx 0.0134          5
# ##  2 dppR02  R02TMCRET mndist mndistRank0 NA      NA           maxdist NA           NA     NA    mnsqrtdist NA    NA    0.00254         5
#
# resdppR01a <- res %>% filter(Outcome=="dppR01") %>% filter(num_stats>=5) %>% filter(p < .03)
# resdppR01long <- pivot_longer(resdppR01a,cols=all_of(outcome_names))
# sort(table(resdppR01long$value,exclude=c()))
# nrow(resdppR01long)
# nrow(resdppR01a)
# resdppR01a %>% summarize(across(one_of(outcome_names),num_not_na))
# #      Y mndist mndistRank0 maddist maddistRank0 maxdist maxdistRank0 mhdist rankx mnsqrtdist hubmn tanhx
# #  <int>  <int>       <int>   <int>        <int>   <int>        <int>  <int> <int>      <int> <int> <int>
# #1    67     22          43      10           13      56           15     21    32         21    15    67
# ## So: Always Y,tanhx
# ## Next best maxdist,rankx,mndistRank0(worst are maddist and maddistRank0 and maxdistRank0 and hubmn)
# ## Intermediate mhdist, mnsqrtdist
#
# resdppR01a %>% filter(!is.na(maxdist) & !is.na(rankx) &  !is.na(mndistRank0))
#
# good_stats_nms <- c("Y","tanhx","mndist","mndistRank0","maxdist","rankx")
#
# resdppR01long %>% filter(value %in% good_stats_nms)
#
# res2 <- res %>% rowwise() %>% mutate(good_stats = num_not_na(c_across(one_of(good_stats_nms)))) %>% ungroup()
# res2 %>% filter(good_stats == length(good_stats_nms) & num_stats==length(good_stats_nms))
#
# ## cands <- res %>% filter(is.na(mndist) & is.na(mndistRank0) & is.na(zscoreY) & is.na(maddistRank0) & is.na(maxdistRank0)) %>%
# ## 	group_by(Outcome) %>%
# ## 	filter(p <= .03) %>%
# ## 	filter(num_stats >= 5) %>%
# ## 	filter(p==min(p)) %>%
# ## 	filter(num_stats==max(num_stats))
# ## cands
# ##
# ##
# ## cands %>% group_by(Outcome) %>% summarize(across(one_of(outcome_names),num_not_na),n=n())
# ##
# ##
# ## dat[,R01md:= (R01TMCRET - mean(R01TMCRET)),by=blockF]
# ##
# ## wilcox_test(R01md~trtF,data=dat)
# ## wilcox_test(R01md~trtF|blockF,data=dat)
# ##
# ## ## Not sure that it matters if we have some blocks with constant outcomes.
# ## dppb1 <- dat[,.(avgR01mdsd=sd(R01md),avgR01sd=sd(R01TMCRET)),by=blockF]
# ## dppb1[avgR01sd==0,]
# ##
# ## dist_lst <- lapply(split(dat,dat$blockF),function(dat){ vecdist(dat$R01TMCRET) })

