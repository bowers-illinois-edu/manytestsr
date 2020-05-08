# Test and develop functions to create treatment effects in simulations
# which blocks did it choose? Did it choose the right blocks?

# context("Simulation Functions")

# library(here)
# source(here::here("tests/testthat", "make_test_data.R"))
# devtools::load_all() ## comment this out for production
setDTthreads(1)
options(digits = 4)
#####
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]
setkey(bdat4, bF)
setkey(bdat, bF)

set.seed(12355)
idat$y1test_null <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 0, prop_blocks_0 = .5)
idat$y1test_zeros <- create_effects(idat = idat, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 2, prop_blocks_0 = .5)
idat[, Y_zeros := y1test_zeros * Z + y0 * (1 - Z)]
idat[, Y_null := y1test_null * Z + y0 * (1 - Z)]

### Add testthat here for these two to assess the create_effects function
## Should not reject
oneway_test(Y_null ~ ZF | bF, data = idat)
## Should reject: half blocks are 0 but they are big enough and the effects are large enough

alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "NULL"),
  sfn = c(
    "splitCluster", "splitEqualApprox", "splitLOO",
    "splitSpecifiedFactor"
  ),
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <- "lvs"
oneway_test(Y_zeros ~ ZF | bF, data = idat)
####

## Now using the function:


err_testing_fn <- function(afn, sfn, sby, fmla = Ytauv2 ~ ZF | bF, idat = idat3, bdat = bdat4, truevar_name) {
  if (afn == "NULL") {
    theafn <- NULL
  } else {
    theafn <- get(afn)
  }
  ## afn and sfn and sby are character names
  theres <- findBlocks(
    idat = idat, bdat = bdat, blockid = "bF", splitfn = get(sfn),
    pfn = pIndepDist, alphafn = theafn, thealpha = 0.05,
    fmla = fmla, # Ynorm_inc ~ ZF | bF,
    parallel = "no", copydts = TRUE, splitby = sby
  )

  detects <- report_detections(theres, only_hits = FALSE, fwer = is.null(afn))

  nodes <- detects[, .(
    hit = as.numeric(unique(hit)),
    numnull = sum(get(truevar_name) == 0),
    anynull = as.numeric(any(get(truevar_name) == 0)),
    allnull = as.numeric(all(get(truevar_name) == 0)),
    numnotnull = sum(get(truevar_name) != 0),
    anynotnull = as.numeric(any(get(truevar_name) != 0))
  ), by = hit_grp]

  # err_tab1 <- with(nodes1 , table(hit, anynotnull, exclude = c()))
  err_tab0 <- with(nodes, table(hit, allnull, exclude = c()))

  if (!identical(dim(err_tab0), as.integer(c(2, 2)))) {
    blank_mat <- matrix(0, 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
    blank_mat[rownames(err_tab0), colnames(err_tab0)] <- err_tab0
    err_tab <- blank_mat
  } else {
    err_tab <- err_tab0
  }

  errs <- calc_errs(theres,
    truevar_name = truevar_name,
    trueeffect_tol = .Machine$double.eps
  )

  expect_equal(errs$true_pos_prop, err_tab["1", "0"] / sum(err_tab))
  expect_equal(errs$true_neg_prop, err_tab["0", "1"] / sum(err_tab))
  expect_equal(errs$false_pos_prop, err_tab["1", "1"] / sum(err_tab))
  expect_equal(errs$false_neg_prop, err_tab["0", "0"] / sum(err_tab))
  expect_equal(errs$true_disc_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))
  expect_equal(errs$false_disc_prop, err_tab["1", "1"] / max(1, sum(err_tab["1", ])))
  expect_equal(errs$true_nondisc_prop, err_tab["0", "1"] / max(1, sum(err_tab["0", ])))
  expect_equal(errs$false_nondisc_prop, err_tab["0", "0"] / max(1, sum(err_tab["0", ])))
}


err_testing_fn2 <- function(fmla = Ytauv2 ~ ZF, idat = idat3, bdat = bdat4, truevar_name, thealpha = .05) {
  ## afn and sfn and sby are character names
  theres <- adjust_block_tests(
    idat = idat, bdat = bdat, blockid = "bF", p_adj_method = "BH",
    pfn = pIndepDist, fmla = fmla, copydts = TRUE
  )

  theres[, hit := max_p < .05] ## doing FDR==.05 for now. All that really matters is some fixed number.

  nodes <- theres[, .(
    bF = bF,
    hit = as.numeric(hit),
    anynull = as.numeric(get(truevar_name) == 0),
    allnull = as.numeric(get(truevar_name) == 0),
    anynotnull = as.numeric(get(truevar_name) != 0)
  )]

  # err_tab1 <- with(nodes1 , table(hit, anynotnull, exclude = c()))
  err_tab0 <- with(nodes, table(hit, allnull, exclude = c()))

  if (!identical(dim(err_tab0), as.integer(c(2, 2)))) {
    blank_mat <- matrix(0, 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
    blank_mat[rownames(err_tab0), colnames(err_tab0)] <- err_tab0
    err_tab <- blank_mat
  } else {
    err_tab <- err_tab0
  }

  errs <- calc_errs(
    testobj = theres,
    truevar_name = truevar_name,
    trueeffect_tol = .Machine$double.eps
  )

  expect_equal(errs$true_pos_prop, err_tab["1", "0"] / sum(err_tab))
  expect_equal(errs$true_neg_prop, err_tab["0", "1"] / sum(err_tab))
  expect_equal(errs$false_pos_prop, err_tab["1", "1"] / sum(err_tab))
  expect_equal(errs$false_neg_prop, err_tab["0", "0"] / sum(err_tab))
  expect_equal(errs$true_disc_prop, err_tab["1", "0"] / max(1, sum(err_tab["1", ])))
  expect_equal(errs$false_disc_prop, err_tab["1", "1"] / max(1, sum(err_tab["1", ])))
  expect_equal(errs$true_nondisc_prop, err_tab["0", "1"] / max(1, sum(err_tab["0", ])))
  expect_equal(errs$false_nondisc_prop, err_tab["0", "0"] / max(1, sum(err_tab["0", ])))
}

# err_testing_fn2(fmla = Ytauv2 ~ ZF, idat = idat3, bdat = bdat4, truevar_name="ate_tauv2")

resnms <- apply(alpha_and_splits, 1, function(x) {
  paste(x, collapse = "_", sep = "")
})


test_that("Error calculations for a given set of tests work: No effects at all", {
  ### No effects at all
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_null") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Ynull ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  ## names(tau_null) <- resnms
})

test_that("Error calculations for a given set of tests work: large and homogenous effects", {
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_homog") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Yhomog ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
})


test_that("Error calculations for a given set of tests work: individually heteogeneous effects and increase with block size. Also some completely null blocks.", {
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_norm_inc") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Ynorm_inc ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
})


test_that("Error calculations for a given set of tests work:individually heteogeneous effects and decrease with block size. Also some completely null blocks.", {
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_norm_dec") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Ynorm_dec ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
})

test_that("Error calculations for a given set of tests work:constant effect that cancel out at the high level (half large and positive, half large and negative).", {
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_tau") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Y ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
})



test_that("Error calculations for a given set of tests work:individually heterogeneous effects block-fixed effects and some completely null blocks.", {
  test_lst <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby, truevar_name = "ate_tauv2") {
      message(paste(afn, sfn, sby, collapse = ","))
      err_testing_fn(
        afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4,
        fmla = Ytauv2 ~ ZF | bF, truevar_name = truevar_name
      )
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
})


## START HERE: Add bottom up testing.
err_testing_fn2(fmla = Ynorm_inc ~ ZF, idat = idat3, bdat = bdat4, truevar_name = "ate_norm_inc")

test_that("Error calculations for a given set of tests work: Testing in every block and adjusting p-values via FDR/BH.", {
  truevar_names <- grep("^ate", names(bdat4), value = TRUE)
  outcome_names <- paste0("Y", gsub("ate_", "", truevar_names))
  outcome_names[1] <- "Y"
  stopifnot(all(outcome_names %in% names(idat3)))

  test_lst <- mapply(
    FUN = function(truevar_name, outcome_name) {
      fmla <- reformulate("ZF", response = outcome_name)
      message(paste(truevar_name, collapse = ","))
      err_testing_fn2(
        idat = idat3, bdat = bdat4,
        fmla = fmla, truevar_name = truevar_name
      )
    },
    truevar_name = truevar_names,
    outcome_name = outcome_names,
    SIMPLIFY = FALSE
  )
})


###
### Test the padj sim function that we use.
## p_sims_res <- parSapplyLB(cl,1:nrow(simparms), FUN=function(i) {
set.seed(12345)
i <- 1
x <- alpha_and_splits[i, ]
xnm <- paste(x, collapse = "_")
message(xnm)
nsims <- 10

p_sims_tab <- padj_test_fn(
  idat = idat3,
  bdat = bdat4,
  blockid = "bF",
  trtid = "Z",
  fmla = Ynorm_inc ~ ZF | blockF,
  ybase = "y0",
  prop_blocks_0 = .5,
  tau_fn = tau_norm_covariate_outliers,
  tau_size = .5,
  covariate = "v4",
  pfn = pIndepDist,
  afn = getFromNamespace(x[["afn"]], ns = "manytestsr"),
  nsims = nsims,
  ncores = 1, ## parallelize over the  rows of simparms
  reveal_and_test_fn = reveal_po_and_test_siup,
  splitfn = getFromNamespace(x[["sfn"]], ns = "manytestsr"),
  splitby = x[["splitby"]]
)

## Testing the error simulation function (deprecated I think. Not used anymore)
## set.seed(12345)
## errsimres1 <- errsimfn(
##   idat = idat, bdat = bdat, afn = NULL, pfn = pIndepDist, splitfn = splitLOO,
##   fmla = Y ~ ZF | bF,
##   block = "bF",
##   trt = "Z", ## must be a binary numeric variable
##   ybase = "y0",
##   thealpha = 0.05, tau_fn = tau_norm, tau_size = 5, prop_blocks_0 = .5,
##   sims = 1000
## )
## set.seed(12345)
## errsimres2 <- errsimfn(
##   idat = idat, bdat = bdat, afn = alpha_saffron, pfn = pIndepDist, splitfn = splitLOO,
##   fmla = Y ~ ZF | bF,
##   block = "bF",
##   trt = "Z", ## must be a binary numeric variable
##   ybase = "y0",
##   thealpha = 0.05, tau_fn = tau_norm, tau_size = 5, prop_blocks_0 = .5,
##   sims = 1000
## )
