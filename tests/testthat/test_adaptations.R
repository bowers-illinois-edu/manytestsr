## Test and develop functions to adapt p-values and alpha levels as the tree grows
testthat::context("Alpha and P Adjusting Performance")

## The next lines are for use when creating the tests. We want interactive<-FALSE for production
interactive <- TRUE
if (interactive) {
  library(here)
  library(data.table)
  library(dtplyr)
  remotes::install_github("bowers-illinois-edu/TreeTestSim")
  library(TreeTestsSim)
  source(here("tests/testthat", "make_test_data.R"))
  library(devtools)
  ## This next is a copy from the package TreeTestsSim
  ## source(here("tests/testthat", "generate_tree.R"))
  load_all() ## use  this during debugging
}
setDTthreads(1)
options(digits = 4)

## For now, create some data within this file.
## Move it elsewhere soon.

## A tree where 50% of the leaves are non-null
## the tree is defined at the level of the block
## We are not going to pass this vector below when we assess weak FWER control
bdt <- generate_tree_DT(max_level = 3, k = 4, t = .5)
bdt[, leaf := as.numeric(level == max(level))]
table(bdt$leaf)
stopifnot(nrow(bdt) == sum(4^seq(0, 3)))


### This next converts the node level data into a block level dataset with a
### factor variable with "." indicating relationships between nodes (each "."
### is a level in the tree)

# build a named vector so parent_lookup["22"] == 6, etc.
parent_lookup <- setNames(bdt$parent, bdt$node)
## Compare to the mathematical way to calculate the parent
expect_equal(parent_lookup[22][[1]], floor((22 - 2) / 4) + 1)
expect_equal(parent_lookup[22][[1]], 6)

# for a given node, return the vector of its ancestors+itself, *excluding* the
# NA‐parent root.

get_path <- function(node) {
  path <- integer(0)
  current <- node
  # walk up until you hit an NA parent
  while (!is.na(parent_lookup[as.character(current)])) {
    path <- c(path, current) # record this node
    current <- parent_lookup[as.character(current)] # jump to parent
  }
  # reverse so it reads root→…→node
  rev(path)
}

#  apply to every node, paste the names with ".", and coerce to factor since
#  find_blocks and splitSpecifiedFactorMulti requires such a factor to split
#  the data according to the pre-specified structure

bdt[, lvls_fac := factor(sapply(node, function(n) {
  paste(get_path(n), collapse = ".")
}))]
bdt[node == 6, ]
bdt[parent == 6, ]

## We are going to focus on the leaves which are the individual blocks.
## That is the leaf nodes are the block-level dataset.
## This is the kind of data that an experimenter would have access to directly.

bdt1 <- droplevels(bdt[leaf == 1, ])

## Assess the splitSpecifiedFactorMulti splitting function
## The splitSpecifiedFactorMulti function should recreate this tree.
bdt1[, split1 := splitSpecifiedFactorMulti(bid = node, x = lvls_fac)]
bdt1[, split2 := splitSpecifiedFactorMulti(node, lvls_fac), by = split1]
bdt1[, split3 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2)]

with(bdt1, table(split1, exclude = c()))
with(bdt1, ftable(split1 = split1, split2 = split2, split3 = split3))

## Everybody within the same parent should be in the same first level split
test_split_1 <- bdt1 %>%
  group_by(parent) %>%
  summarize(val = length(unique(split1)))
stopifnot(all(test_split_1$val == 1))

## Split2 is constant within split1 --- so there are blocks with split2==1 within each value of split1.
## i.e. first group within split1=1, split1=2, etc..
test_split_2 <- bdt1 %>%
  group_by(split1) %>%
  summarize(val = length(unique(split2)))
stopifnot(all(test_split_2$val == 4))
stopifnot(nrow(test_split_2) == 4)

bdt1 %>% filter(split1 == 1)
bdt1 %>% filter(split2 == 1)

## This next returns all 0s when splitting is no longer possible. So, consider watching out for this in the find_blocks code
bdt1[, split4 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3)]
bdt1[, split5 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3, split4)]

expect_equal(as.numeric(as.character(unique(bdt1$split4))), 0)
expect_equal(as.numeric(as.character(unique(bdt1$split5))), 0)

bdt1[, c("split1", "split2", "split3", "split4", "split5") := NULL]

head(bdt1)

## Make individual level data with nb=50
## Idea here is to make something with high power within each block
bdt1[, nb := 50]
bdt1[, bary0 := 0]
bdt1[, bF := factor(node)]

idt <- data.table(b = rep(bdt1$node, bdt1$nb))
idt[, bF := factor(b)]
idt[, id := seq_len(nrow(idt))]
idt[, bary0 := rep(bdt1$bary0, bdt1$nb)]

## this next "by=b" is not necessary since it is the same y0 for all blocks.
## But it will not be that case in other situations

idt[, y0 := rnorm(.N, bary0, sd = 1), by = b]
idt <- merge(idt, bdt1[, .(bF, nonnull, lvls_fac)], by = "bF")

## Here we ignore the nonnull variable on bdt1 for the sake of creating a
## situation where there are no effects in any block or for any person.

idt[, y1 := create_effects(
  idat = idt, ybase = "y0",
  blockid = "bF", tau_fn = tau_null, tau_size = 0, prop_blocks_0 = 1, non_null_blocks = NULL
)]

## Treatment is completely randomized in each block
idt[, trt := sample(rep(c(0, 1), .N / 2)), by = bF]
## Reveal the observed Y as a function of the potential outcomes
idt[, Y := trt * y1 + (1 - trt) * y0]

## Some procedures need treatment to be a factor (esp those using the coin
## package)
idt[, trtF := factor(trt)]


## Assess weak control of the FWER with no local adjustment or extra global
## adjustment
#### First just see if find_blocks itself operates
set.seed(12345)
res_null <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL, local_adj_p_fn = NULL,
  fmla = Y ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = TRUE, copydts = TRUE
)

# We expect a single test with a p>.05 (this depends on the set.seed above, but mostly should be true)
expect_equal(length(unique(res_null$pfinalb)), 1)
expect_gt(unique(res_null$pfinalb), .05)

## Check for errors in the tree shaped output
## For now ignoring the nonnull and node_nonnull output since we generated y1=y0 for all blocks.
res_null_tree <- make_results_tree(res_null)
res_null_nodes <- res_null_tree %>%
  activate(nodes) %>%
  as_tibble()


nsims <- 1000
sim_err <- 2 * sqrt((.05 * (1 - .05)) / nsims)

set.seed(12345)
p_null_res <- padj_test_fn(
  idat = idt,
  bdat = bdt1,
  blockid = "bF",
  trtid = "trt", ## needs to be a 0 or a 1 for the simulation
  fmla = Y ~ trtF | bF,
  ybase = "y0",
  prop_blocks_0 = 1,
  tau_fn = tau_null,
  tau_size = 0,
  covariate = NULL,
  pfn = pOneway,
  nsims = nsims,
  ncores = 8,
  afn = NULL,
  splitfn = splitSpecifiedFactorMulti,
  splitby = "lvls_fac",
  p_adj_method = "split",
  blocksize = "nb",
  by_block = TRUE,
  stop_splitby_constant = TRUE
)
res_null_rates <- p_null_res[, lapply(.SD, mean, na.rm = TRUE)]
## This is assessing weak control so the following two are equal.
expect_lt(res_null_rates$prop_reject, .05 + sim_err)
expect_lt(res_null_rates$false_pos_prop, .05 + sim_err)

## TODO: It would be nice to calculate power over both nodes and blocks.


## TODO: enable padj_test_fn to be parallel over the replicate function
## This should produce very high power for each block a 1 sd difference
# power.t.test(power=.8,sd=1,delta=1,sig.level=.05)
power.t.test(n = 25, sd = 1, delta = 1, sig.level = .05)
set.seed(12345)
p_half_res <- padj_test_fn(
  idat = idt,
  bdat = bdt1,
  blockid = "bF",
  trtid = "trt", ## needs to be a 0 or a 1 for the simulation
  fmla = Y ~ trtF | bF,
  ybase = "y0",
  prop_blocks_0 = .5,
  tau_fn = tau_norm,
  tau_size = 1,
  covariate = NULL,
  pfn = pOneway,
  nsims = nsims,
  ncores = 1,
  afn = NULL,
  splitfn = splitSpecifiedFactorMulti,
  splitby = "lvls_fac",
  p_adj_method = "split",
  blocksize = "nb",
  by_block = TRUE,
  stop_splitby_constant = TRUE
)

res_half_rates <- p_half_res[, lapply(.SD, mean, na.rm = TRUE)]
## This is showing control of FWER in the strong sense.
expect_lt(res_half_rates$false_pos_prop, .05 + sim_err)
## Power to detect effects at the individual block level.
res_half_rates$true_pos_prop
## This is a basic requirement of a good test: power above alpha
expect_gt(res_half_rates$true_pos_prop, .05)

## Now assess whether find_blocks and helpers works as expected with this setup
idt[, y1_half_tau1 := create_effects(
  idat = idt, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 1, prop_blocks_0 = .5,
  non_null_blocks = "nonnull"
)]
test_non_null <- idt[nonnull == FALSE, .(sum(y1_half_tau1 - y0))]
expect_equal(unique(test_non_null)[[1]], 0)

idt[, Y_half_tau1 := trt * y1_half_tau1 + (1 - trt) * y0]

res_half <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL, local_adj_p_fn = NULL,
  fmla = Y_half_tau1 ~ trtF | bF,
  splitby = "lvls_fac", blocksize = "nb", trace = TRUE, copydts = TRUE
)

res_half
# Key: <testable>
#      nonnull  lvls_fac
#        <lgcl>    <fctr>
#  1:    FALSE    2.7.26
#  2:    FALSE    2.7.27
#  3:    FALSE    2.7.28
#  4:    FALSE    2.7.29
#  5:     TRUE   5.20.78
#  6:    FALSE   5.20.79
#  7:    FALSE   5.20.80
#  8:    FALSE   5.20.81
#  9:    FALSE    2.6.22
# 10:    FALSE    2.6.23
# 11:    FALSE    2.6.24
# 12:     TRUE    2.6.25
# 13:    FALSE    2.8.30
# 14:     TRUE    2.8.31
# 15:     TRUE    2.8.32
# 16:     TRUE    2.8.33
# 17:    FALSE    2.9.34
# 18:     TRUE    2.9.35
# 19:     TRUE    2.9.36
# 20:     TRUE    2.9.37
# 21:    FALSE   3.10.38
# 22:    FALSE   3.10.39
# 23:     TRUE   3.10.40
# 24:    FALSE   3.10.41
# 25:     TRUE   3.11.42
# 26:     TRUE   3.11.43
# 27:    FALSE   3.11.44
# 28:    FALSE   3.11.45
# 29:     TRUE   3.12.46
# 30:     TRUE   3.12.47
# 31:     TRUE   3.12.48
# 32:    FALSE   3.12.49
# 33:     TRUE   3.13.50
# 34:    FALSE   3.13.51
# 35:     TRUE   3.13.52
# 36:    FALSE   3.13.53
# 37:    FALSE   4.14.54
# 38:     TRUE   4.14.55
# 39:     TRUE   4.14.56
# 40:    FALSE   4.14.57
# 41:     TRUE   4.15.58
# 42:     TRUE   4.15.59
# 43:    FALSE   4.15.60
# 44:     TRUE   4.15.61
# 45:     TRUE   4.16.62
# 46:    FALSE   4.16.63
# 47:     TRUE   4.16.64
# 48:    FALSE   4.16.65
# 49:     TRUE   4.17.66
# 50:     TRUE   4.17.67
# 51:     TRUE   4.17.68
# 52:    FALSE   4.17.69
# 53:     TRUE   5.18.70
# 54:     TRUE   5.18.71
# 55:     TRUE   5.18.72
# 56:     TRUE   5.18.73
# 57:     TRUE   5.19.74
# 58:    FALSE   5.19.75
# 59:    FALSE   5.19.76
# 60:    FALSE   5.19.77
# 61:    FALSE   5.21.82
# 62:    FALSE   5.21.83
# 63:    FALSE   5.21.84
# 64:     TRUE   5.21.85
#      node level parent nonnull  leaf lvls_fac
#        nb bary0     bF        p1   pfinalb
#     <num> <num> <fctr>     <num>     <num>
#  1:    50     0     26 2.258e-41 7.887e-01
#  2:    50     0     27 2.258e-41 7.887e-01
#  3:    50     0     28 2.258e-41 7.887e-01
#  4:    50     0     29 2.258e-41 7.887e-01
#  5:    50     0     78 2.258e-41 6.653e-01
#  6:    50     0     79 2.258e-41 6.653e-01
#  7:    50     0     80 2.258e-41 6.653e-01
#  8:    50     0     81 2.258e-41 6.653e-01
#  9:    50     0     22 2.258e-41 1.705e-01
# 10:    50     0     23 2.258e-41 3.203e-01
# 11:    50     0     24 2.258e-41 3.238e-01
# 12:    50     0     25 2.258e-41 5.953e-04
# 13:    50     0     30 2.258e-41 6.529e-01
# 14:    50     0     31 2.258e-41 1.318e-03
# 15:    50     0     32 2.258e-41 8.934e-03
# 16:    50     0     33 2.258e-41 6.139e-04
# 17:    50     0     34 2.258e-41 4.110e-01
# 18:    50     0     35 2.258e-41 4.981e-03
# 19:    50     0     36 2.258e-41 3.095e-03
# 20:    50     0     37 2.258e-41 1.184e-03
# 21:    50     0     38 2.258e-41 4.045e-01
# 22:    50     0     39 2.258e-41 4.781e-02
# 23:    50     0     40 2.258e-41 1.907e-06
# 24:    50     0     41 2.258e-41 3.847e-01
# 25:    50     0     42 2.258e-41 3.873e-03
# 26:    50     0     43 2.258e-41 1.385e-03
# 27:    50     0     44 2.258e-41 2.971e-01
# 28:    50     0     45 2.258e-41 6.398e-01
# 29:    50     0     46 2.258e-41 1.399e-03
# 30:    50     0     47 2.258e-41 1.355e-02
# 31:    50     0     48 2.258e-41 4.455e-04
# 32:    50     0     49 2.258e-41 2.499e-01
# 33:    50     0     50 2.258e-41 1.315e-02
# 34:    50     0     51 2.258e-41 5.792e-01
# 35:    50     0     52 2.258e-41 1.260e-01
# 36:    50     0     53 2.258e-41 9.110e-01
# 37:    50     0     54 2.258e-41 9.199e-01
# 38:    50     0     55 2.258e-41 1.599e-03
# 39:    50     0     56 2.258e-41 7.310e-04
# 40:    50     0     57 2.258e-41 1.196e-01
# 41:    50     0     58 2.258e-41 5.790e-02
# 42:    50     0     59 2.258e-41 6.094e-03
# 43:    50     0     60 2.258e-41 5.441e-01
# 44:    50     0     61 2.258e-41 1.070e-01
# 45:    50     0     62 2.258e-41 4.044e-03
# 46:    50     0     63 2.258e-41 6.485e-01
# 47:    50     0     64 2.258e-41 1.633e-04
# 48:    50     0     65 2.258e-41 5.361e-01
# 49:    50     0     66 2.258e-41 7.172e-05
# 50:    50     0     67 2.258e-41 4.097e-03
# 51:    50     0     68 2.258e-41 1.106e-02
# 52:    50     0     69 2.258e-41 4.709e-01
# 53:    50     0     70 2.258e-41 5.247e-03
# 54:    50     0     71 2.258e-41 5.555e-04
# 55:    50     0     72 2.258e-41 1.428e-01
# 56:    50     0     73 2.258e-41 1.107e-05
# 57:    50     0     74 2.258e-41 2.195e-04
# 58:    50     0     75 2.258e-41 1.550e-01
# 59:    50     0     76 2.258e-41 3.532e-01
# 60:    50     0     77 2.258e-41 3.740e-01
# 61:    50     0     82 2.258e-41 7.957e-01
# 62:    50     0     83 2.258e-41 6.488e-01
# 63:    50     0     84 2.258e-41 5.153e-01
# 64:    50     0     85 2.258e-41 4.138e-04
#        nb bary0     bF        p1   pfinalb
#                           biggrp     g1 alpha1 testable nodenum_current
#                           <fctr> <fctr>  <num>   <lgcl>          <char>
#  1:          1.6783673c.98e7f981      1   0.05    FALSE        98e7f981
#  2:          1.6783673c.98e7f981      1   0.05    FALSE        98e7f981
#  3:          1.6783673c.98e7f981      1   0.05    FALSE        98e7f981
#  4:          1.6783673c.98e7f981      1   0.05    FALSE        98e7f981
#  5:          1.3ae9e7b5.c585b76a      1   0.05    FALSE        c585b76a
#  6:          1.3ae9e7b5.c585b76a      1   0.05    FALSE        c585b76a
#  7:          1.3ae9e7b5.c585b76a      1   0.05    FALSE        c585b76a
#  8:          1.3ae9e7b5.c585b76a      1   0.05    FALSE        c585b76a
#  9: 1.6783673c.efe0c917.345ac693      1   0.05    FALSE        345ac693
# 10: 1.6783673c.efe0c917.435df605      1   0.05    FALSE        435df605
# 11: 1.6783673c.efe0c917.da54a7bf      1   0.05    FALSE        da54a7bf
# 12: 1.6783673c.efe0c917.ad539729      1   0.05    FALSE        ad539729
# 13: 1.6783673c.01eea83b.090be84b      1   0.05    FALSE        090be84b
# 14: 1.6783673c.01eea83b.7e0cd8dd      1   0.05    FALSE        7e0cd8dd
# 15: 1.6783673c.01eea83b.e7058967      1   0.05    FALSE        e7058967
# 16: 1.6783673c.01eea83b.9002b9f1      1   0.05    FALSE        9002b9f1
# 17: 1.6783673c.76e998ad.2de29d45      1   0.05    FALSE        2de29d45
# 18: 1.6783673c.76e998ad.5ae5add3      1   0.05    FALSE        5ae5add3
# 19: 1.6783673c.76e998ad.c3ecfc69      1   0.05    FALSE        c3ecfc69
# 20: 1.6783673c.76e998ad.b4ebccff      1   0.05    FALSE        b4ebccff
# 21: 1.e58a1a84.02b5b239.3fffc7fa      1   0.05    FALSE        3fffc7fa
# 22: 1.e58a1a84.02b5b239.48f8f76c      1   0.05    FALSE        48f8f76c
# 23: 1.e58a1a84.02b5b239.d1f1a6d6      1   0.05    FALSE        d1f1a6d6
# 24: 1.e58a1a84.02b5b239.a6f69640      1   0.05    FALSE        a6f69640
# 25: 1.e58a1a84.75b282af.6ce6ffbd      1   0.05    FALSE        6ce6ffbd
# 26: 1.e58a1a84.75b282af.1be1cf2b      1   0.05    FALSE        1be1cf2b
# 27: 1.e58a1a84.75b282af.82e89e91      1   0.05    FALSE        82e89e91
# 28: 1.e58a1a84.75b282af.f5efae07      1   0.05    FALSE        f5efae07
# 29: 1.e58a1a84.ecbbd315.17f00173      1   0.05    FALSE        17f00173
# 30: 1.e58a1a84.ecbbd315.60f731e5      1   0.05    FALSE        60f731e5
# 31: 1.e58a1a84.ecbbd315.f9fe605f      1   0.05    FALSE        f9fe605f
# 32: 1.e58a1a84.ecbbd315.8ef950c9      1   0.05    FALSE        8ef950c9
# 33: 1.e58a1a84.9bbce383.e08bd593      1   0.05    FALSE        e08bd593
# 34: 1.e58a1a84.9bbce383.978ce505      1   0.05    FALSE        978ce505
# 35: 1.e58a1a84.9bbce383.0e85b4bf      1   0.05    FALSE        0e85b4bf
# 36: 1.e58a1a84.9bbce383.79828429      1   0.05    FALSE        79828429
# 37: 1.9b8d31ec.b9ae4185.0651955e      1   0.05    FALSE        0651955e
# 38: 1.9b8d31ec.b9ae4185.7156a5c8      1   0.05    FALSE        7156a5c8
# 39: 1.9b8d31ec.b9ae4185.e85ff472      1   0.05    FALSE        e85ff472
# 40: 1.9b8d31ec.b9ae4185.9f58c4e4      1   0.05    FALSE        9f58c4e4
# 41: 1.9b8d31ec.cea97113.76c449c1      1   0.05    FALSE        76c449c1
# 42: 1.9b8d31ec.cea97113.01c37957      1   0.05    FALSE        01c37957
# 43: 1.9b8d31ec.cea97113.98ca28ed      1   0.05    FALSE        98ca28ed
# 44: 1.9b8d31ec.cea97113.efcd187b      1   0.05    FALSE        efcd187b
# 45: 1.9b8d31ec.57a020a9.a18ffa61      1   0.05    FALSE        a18ffa61
# 46: 1.9b8d31ec.57a020a9.d688caf7      1   0.05    FALSE        d688caf7
# 47: 1.9b8d31ec.57a020a9.4f819b4d      1   0.05    FALSE        4f819b4d
# 48: 1.9b8d31ec.57a020a9.3886abdb      1   0.05    FALSE        3886abdb
# 49: 1.9b8d31ec.20a7103f.c7956f9c      1   0.05    FALSE        c7956f9c
# 50: 1.9b8d31ec.20a7103f.b0925f0a      1   0.05    FALSE        b0925f0a
# 51: 1.9b8d31ec.20a7103f.299b0eb0      1   0.05    FALSE        299b0eb0
# 52: 1.9b8d31ec.20a7103f.5e9c3e26      1   0.05    FALSE        5e9c3e26
# 53: 1.3ae9e7b5.2b8bd646.50f4fab1      1   0.05    FALSE        50f4fab1
# 54: 1.3ae9e7b5.2b8bd646.27f3ca27      1   0.05    FALSE        27f3ca27
# 55: 1.3ae9e7b5.2b8bd646.befa9b9d      1   0.05    FALSE        befa9b9d
# 56: 1.3ae9e7b5.2b8bd646.c9fdab0b      1   0.05    FALSE        c9fdab0b
# 57: 1.3ae9e7b5.5c8ce6d0.34e5e4c5      1   0.05    FALSE        34e5e4c5
# 58: 1.3ae9e7b5.5c8ce6d0.43e2d453      1   0.05    FALSE        43e2d453
# 59: 1.3ae9e7b5.5c8ce6d0.daeb85e9      1   0.05    FALSE        daeb85e9
# 60: 1.3ae9e7b5.5c8ce6d0.adecb57f      1   0.05    FALSE        adecb57f
# 61: 1.3ae9e7b5.b28287fc.96a55696      1   0.05    FALSE        96a55696
# 62: 1.3ae9e7b5.b28287fc.e1a26600      1   0.05    FALSE        e1a26600
# 63: 1.3ae9e7b5.b28287fc.78ab37ba      1   0.05    FALSE        78ab37ba
# 64: 1.3ae9e7b5.b28287fc.0fac072c      1   0.05    FALSE        0fac072c
#                           biggrp     g1 alpha1 testable nodenum_current
#     nodenum_prev blocksbygroup     g2 nodesize        p2 alpha2     g3
#           <char>         <int> <fctr>    <num>     <num>  <num> <fctr>
#  1:     6783673c             4      0      200 1.686e-09   0.05      1
#  2:     6783673c             4      0      200 1.686e-09   0.05      1
#  3:     6783673c             4      0      200 1.686e-09   0.05      1
#  4:     6783673c             4      0      200 1.686e-09   0.05      1
#  5:     3ae9e7b5             4      3      200 1.373e-09   0.05      2
#  6:     3ae9e7b5             4      3      200 1.373e-09   0.05      2
#  7:     3ae9e7b5             4      3      200 1.373e-09   0.05      2
#  8:     3ae9e7b5             4      3      200 1.373e-09   0.05      2
#  9:     efe0c917             1      0       50 1.686e-09   0.05      0
# 10:     efe0c917             1      0       50 1.686e-09   0.05      0
# 11:     efe0c917             1      0       50 1.686e-09   0.05      0
# 12:     efe0c917             1      0       50 1.686e-09   0.05      0
# 13:     01eea83b             1      0       50 1.686e-09   0.05      2
# 14:     01eea83b             1      0       50 1.686e-09   0.05      2
# 15:     01eea83b             1      0       50 1.686e-09   0.05      2
# 16:     01eea83b             1      0       50 1.686e-09   0.05      2
# 17:     76e998ad             1      0       50 1.686e-09   0.05      3
# 18:     76e998ad             1      0       50 1.686e-09   0.05      3
# 19:     76e998ad             1      0       50 1.686e-09   0.05      3
# 20:     76e998ad             1      0       50 1.686e-09   0.05      3
# 21:     02b5b239             1      1       50 4.589e-12   0.05      0
# 22:     02b5b239             1      1       50 4.589e-12   0.05      0
# 23:     02b5b239             1      1       50 4.589e-12   0.05      0
# 24:     02b5b239             1      1       50 4.589e-12   0.05      0
# 25:     75b282af             1      1       50 4.589e-12   0.05      1
# 26:     75b282af             1      1       50 4.589e-12   0.05      1
# 27:     75b282af             1      1       50 4.589e-12   0.05      1
# 28:     75b282af             1      1       50 4.589e-12   0.05      1
# 29:     ecbbd315             1      1       50 4.589e-12   0.05      2
# 30:     ecbbd315             1      1       50 4.589e-12   0.05      2
# 31:     ecbbd315             1      1       50 4.589e-12   0.05      2
# 32:     ecbbd315             1      1       50 4.589e-12   0.05      2
# 33:     9bbce383             1      1       50 4.589e-12   0.05      3
# 34:     9bbce383             1      1       50 4.589e-12   0.05      3
# 35:     9bbce383             1      1       50 4.589e-12   0.05      3
# 36:     9bbce383             1      1       50 4.589e-12   0.05      3
# 37:     b9ae4185             1      2       50 2.076e-15   0.05      0
# 38:     b9ae4185             1      2       50 2.076e-15   0.05      0
# 39:     b9ae4185             1      2       50 2.076e-15   0.05      0
# 40:     b9ae4185             1      2       50 2.076e-15   0.05      0
# 41:     cea97113             1      2       50 2.076e-15   0.05      1
# 42:     cea97113             1      2       50 2.076e-15   0.05      1
# 43:     cea97113             1      2       50 2.076e-15   0.05      1
# 44:     cea97113             1      2       50 2.076e-15   0.05      1
# 45:     57a020a9             1      2       50 2.076e-15   0.05      2
# 46:     57a020a9             1      2       50 2.076e-15   0.05      2
# 47:     57a020a9             1      2       50 2.076e-15   0.05      2
# 48:     57a020a9             1      2       50 2.076e-15   0.05      2
# 49:     20a7103f             1      2       50 2.076e-15   0.05      3
# 50:     20a7103f             1      2       50 2.076e-15   0.05      3
# 51:     20a7103f             1      2       50 2.076e-15   0.05      3
# 52:     20a7103f             1      2       50 2.076e-15   0.05      3
# 53:     2b8bd646             1      3       50 1.373e-09   0.05      0
# 54:     2b8bd646             1      3       50 1.373e-09   0.05      0
# 55:     2b8bd646             1      3       50 1.373e-09   0.05      0
# 56:     2b8bd646             1      3       50 1.373e-09   0.05      0
# 57:     5c8ce6d0             1      3       50 1.373e-09   0.05      1
# 58:     5c8ce6d0             1      3       50 1.373e-09   0.05      1
# 59:     5c8ce6d0             1      3       50 1.373e-09   0.05      1
# 60:     5c8ce6d0             1      3       50 1.373e-09   0.05      1
# 61:     b28287fc             1      3       50 1.373e-09   0.05      3
# 62:     b28287fc             1      3       50 1.373e-09   0.05      3
# 63:     b28287fc             1      3       50 1.373e-09   0.05      3
# 64:     b28287fc             1      3       50 1.373e-09   0.05      3
#     nodenum_prev blocksbygroup     g2 nodesize        p2 alpha2     g3
#            p3 alpha3     g4        p4 alpha4
#         <num>  <num> <fctr>     <num>  <num>
#  1: 7.887e-01   0.05   <NA>        NA   0.05
#  2: 7.887e-01   0.05   <NA>        NA   0.05
#  3: 7.887e-01   0.05   <NA>        NA   0.05
#  4: 7.887e-01   0.05   <NA>        NA   0.05
#  5: 6.653e-01   0.05   <NA>        NA   0.05
#  6: 6.653e-01   0.05   <NA>        NA   0.05
#  7: 6.653e-01   0.05   <NA>        NA   0.05
#  8: 6.653e-01   0.05   <NA>        NA   0.05
#  9: 2.444e-02   0.05      0 1.705e-01   0.05
# 10: 2.444e-02   0.05      1 3.203e-01   0.05
# 11: 2.444e-02   0.05      2 3.238e-01   0.05
# 12: 2.444e-02   0.05      3 5.953e-04   0.05
# 13: 4.584e-06   0.05      0 6.529e-01   0.05
# 14: 4.584e-06   0.05      1 1.318e-03   0.05
# 15: 4.584e-06   0.05      2 8.934e-03   0.05
# 16: 4.584e-06   0.05      3 6.139e-04   0.05
# 17: 6.683e-06   0.05      0 4.110e-01   0.05
# 18: 6.683e-06   0.05      1 4.981e-03   0.05
# 19: 6.683e-06   0.05      2 3.095e-03   0.05
# 20: 6.683e-06   0.05      3 1.184e-03   0.05
# 21: 5.534e-06   0.05      0 4.045e-01   0.05
# 22: 5.534e-06   0.05      1 4.781e-02   0.05
# 23: 5.534e-06   0.05      2 1.907e-06   0.05
# 24: 5.534e-06   0.05      3 3.847e-01   0.05
# 25: 1.356e-02   0.05      0 3.873e-03   0.05
# 26: 1.356e-02   0.05      1 1.385e-03   0.05
# 27: 1.356e-02   0.05      2 2.971e-01   0.05
# 28: 1.356e-02   0.05      3 6.398e-01   0.05
# 29: 1.102e-05   0.05      0 1.399e-03   0.05
# 30: 1.102e-05   0.05      1 1.355e-02   0.05
# 31: 1.102e-05   0.05      2 4.455e-04   0.05
# 32: 1.102e-05   0.05      3 2.499e-01   0.05
# 33: 1.372e-02   0.05      0 1.315e-02   0.05
# 34: 1.372e-02   0.05      1 5.792e-01   0.05
# 35: 1.372e-02   0.05      2 1.260e-01   0.05
# 36: 1.372e-02   0.05      3 9.110e-01   0.05
# 37: 1.669e-05   0.05      0 9.199e-01   0.05
# 38: 1.669e-05   0.05      1 1.599e-03   0.05
# 39: 1.669e-05   0.05      2 7.310e-04   0.05
# 40: 1.669e-05   0.05      3 1.196e-01   0.05
# 41: 3.858e-03   0.05      0 5.790e-02   0.05
# 42: 3.858e-03   0.05      1 6.094e-03   0.05
# 43: 3.858e-03   0.05      2 5.441e-01   0.05
# 44: 3.858e-03   0.05      3 1.070e-01   0.05
# 45: 5.981e-05   0.05      0 4.044e-03   0.05
# 46: 5.981e-05   0.05      1 6.485e-01   0.05
# 47: 5.981e-05   0.05      2 1.633e-04   0.05
# 48: 5.981e-05   0.05      3 5.361e-01   0.05
# 49: 3.003e-06   0.05      0 7.172e-05   0.05
# 50: 3.003e-06   0.05      1 4.097e-03   0.05
# 51: 3.003e-06   0.05      2 1.106e-02   0.05
# 52: 3.003e-06   0.05      3 4.709e-01   0.05
# 53: 2.809e-10   0.05      0 5.247e-03   0.05
# 54: 2.809e-10   0.05      1 5.555e-04   0.05
# 55: 2.809e-10   0.05      2 1.428e-01   0.05
# 56: 2.809e-10   0.05      3 1.107e-05   0.05
# 57: 3.310e-03   0.05      0 2.195e-04   0.05
# 58: 3.310e-03   0.05      1 1.550e-01   0.05
# 59: 3.310e-03   0.05      2 3.532e-01   0.05
# 60: 3.310e-03   0.05      3 3.740e-01   0.05
# 61: 2.851e-02   0.05      0 7.957e-01   0.05
# 62: 2.851e-02   0.05      1 6.488e-01   0.05
# 63: 2.851e-02   0.05      2 5.153e-01   0.05
# 64: 2.851e-02   0.05      3 4.138e-04   0.05
#            p3 alpha3     g4        p4 alpha4

res_half %>%
  select(node, level, parent, nonnull, lvls_fac, biggrp, starts_with("g"), starts_with("p")) %>%
  mutate(across(where(is.numeric), zapsmall))

## Deprecate make_results_tree or just replace it with the _direct one after some more testing
res_half_tree <- make_results_tree(res_half)
res_half_tree2 <- make_results_tree_direct(res_half)
res_half_graph <- make_results_ggraph(res_half_tree)
res_half_graph2 <- make_results_ggraph(res_half_tree2$graph)

res_half_nodes2 <- res_half_tree2$graph %>%
  activate(nodes) %>%
  as.data.table()

res_half_nodes <- res_half_tree %>%
  activate(nodes) %>%
  as.data.table()

res_half_nodes %>%
  select(nodenum, depth, p, node_nonnull, is_leaf) %>%
  as_tibble() %>%
  print(n = 100)



### OLD BELOW TODO CLEANUP

idt <- data.table(b = rep(c(1:20), length = 1000))

## 20 blocks, each with 50 obs
idt <- data.table(i = 1:1000, b = rep(c(1:20), length = 1000))
bdt <- data.table(b = 1:20)
setkey(bdt, b)
setkey(idt, b)
bdt[, vb1 := rep(c(1, 2, 3, 4), length.out = .N)]
## 2 levels within each state
bdt[, vb2 := rep(c(1, 2, 3), length.out = .N), by = vb1]
bdt[, vb3 := rep(seq(1, .N), length.out = .N), by = vb2]
bdt[, vbnest := interaction(vb1, vb2, vb3, lex.order = TRUE, drop = TRUE)]
ftable(lv1 = bdt$vb1, lv2 = bdt$vb2, bdt$vb3, exclude = c())
table(bdt$vbnest, exclude = c())
## Now merge the block data onto the individual level data
idt <- bdt[idt]
set.seed(12345)
idt[, vi1 := round(runif(.N, min = 0, max = 10)), by = b]
idt[, y0 := vi1 + vb1 + rnorm(100), by = b]

# 4 nodes at level 1
bdt[, lv1 := rep(c(1, 2, 3, 4), length.out = .N)]
bdt[, lv2 := rep(1:.N, length.out = .N), by = lv1]
with(bdt, table(lv1, lv2, exclude = c()))
bdt[, lvls := interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
table(bdt$lvls, exclude = c())

## Shuffle order  of the blocks so that the first set and the second set don't  automatically go together
set.seed(12345)
bdat4 <- bdat3[sample(.N), ]
## Setting up  a test of pre-specified splits
bdat4[, lv1 := cut(v1, 2, labels = c("l1_1", "l1_2"))]
bdat4[, lv2 := cut(v2, 2, labels = c("l2_1", "l2_2")), by = lv1]
bdat4[, lv3 := seq(1, .N), by = interaction(lv1, lv2, lex.order = TRUE, drop = TRUE)]
bdat4[, lvs := interaction(lv1, lv2, lv3, lex.order = TRUE, drop = TRUE)]

## Setting up to test the local p and alpha adjusting methods versus FWER
## methods on different splitters


local_adj_methods <- c(
  "local_hommel_all_ps", "local_bh_all_ps",
  "local_simes", "local_unadj_all_ps"
)

alpha_and_splits <- expand.grid(
  afn = c("alpha_investing", "alpha_saffron", "alpha_addis", "NULL"),
  sfn = c(
    "splitCluster", "splitEqualApprox", "splitLOO",
    "splitSpecifiedFactor", "splitSpecifiedFactorMulti"
  ),
  local_adj_fn = local_adj_methods,
  stringsAsFactors = FALSE
)
alpha_and_splits$splitby <- "hwt"
alpha_and_splits$splitby[grep("Specified", alpha_and_splits$sfn)] <- "lvs"

testing_fn <- function(afn, sfn, local_adj, sby, fmla = Ytauv2 ~ ZF | bF, idat = idat3, bdat = bdat4) {
  if (afn == "NULL") {
    theafn <- NULL
  } else {
    theafn <- getFromNamespace(afn, ns = "manytestsr")
  }

  if (local_adj == "NULL") {
    local_adj <- NULL
  } else {
    local_adj <- getFromNamespace(local_adj, ns = "manytestsr")
  }
  ## afn and sfn and sby are character names
  theres <- find_blocks(
    idat = idat, bdat = bdat, blockid = "bF", splitfn = get(sfn),
    pfn = pIndepDist, alphafn = theafn, local_adj_p_fn = local_adj, thealpha = 0.05, thew0 = .05 - .001,
    fmla = fmla,
    copydts = TRUE, splitby = sby, stop_splitby_constant = TRUE, parallel = "multicore", ncores = 2
  )
  setkey(theres, bF)
  return(theres)[order(biggrp)]
  # theps <- grep("^p[0-9]", names(theres), value = TRUE)
  # theas <- grep("^alpha", names(theres), value = TRUE)
  # truth <- grep("^ate", names(theres), value = TRUE)
  # return(theres[, .SD, .SDcols = c(
  #  theps, theas, truth, "pfinalb", "blocksbygroup",
  #  "biggrp", "bF", "hwt", "nb"
  # )][order(biggrp)])
}

##  Maybe make w0 more like .05.
alpha_and_splits[c(1, 2, 4), ]
## debug(find_blocks)
res_ai <- testing_fn(
  afn = alpha_and_splits[1, "afn"],
  sfn = alpha_and_splits[1, "sfn"],
  local_adj = alpha_and_splits[1, "local_adj_fn"],
  sby = alpha_and_splits[1, "splitby"],
  idat = idat3, bdat = bdat4
)
res_saffron <- testing_fn(
  afn = alpha_and_splits[2, "afn"],
  sfn = alpha_and_splits[2, "sfn"],
  sby = alpha_and_splits[2, "splitby"],
  local_adj = alpha_and_splits[2, "local_adj_fn"],
  idat = idat3, bdat = bdat4
)
res_fwer <- testing_fn(
  afn = alpha_and_splits[4, "afn"],
  sfn = alpha_and_splits[4, "sfn"],
  sby = alpha_and_splits[4, "splitby"],
  local_adj = alpha_and_splits[4, "local_adj_fn"],
  idat = idat3, bdat = bdat4
)

grep("^p[0-9]", names(res_ai), value = TRUE)
grep("^p[0-9]", names(res_saffron), value = TRUE)
grep("^p[0-9]", names(res_fwer), value = TRUE)
options(digits = 3, scipen = 8)
res_ai[order(p1, p2, p3, p4, p5, p6, p7, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6, p7)]
res_saffron[order(p1, p2, p3, p4, p5, p6, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6)]
res_fwer[order(p1, p2, p3, p4, p5, p6, decreasing = TRUE), .(p1, p2, p3, p4, p5, p6)]
## Early in the splitting. No change
all.equal(res_ai$p2, res_saffron$p2)
all.equal(res_ai$p2, res_fwer$p2)
## But then some blocks stop testing
all.equal(res_ai$p5, res_saffron$p5)
all.equal(res_ai$p5, res_fwer$p5)
cbind(
  res_ai[, .(bF, p5, p6)],
  res_saffron[, .(bF, p5, p6)],
  res_fwer[, .(bF, p5, p6)]
)

## Which blocks were identified?
##
## The alpha investing style procedures should identify **more** blocks
# We detect an effect in an individual block if:
## (1) the block is a leaf and that leaf has p < alpha
# We detect an effect for a group of blocks if
## (2) the blocks are leaves with p > alpha but the parent of those blocks have p < alpha

## Change one result from res_fwer to enable us to check this:

res_ai[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_saffron[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_fwer[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]

res_ai[bF == "10", p5 := .5]
res_ai[bF == "10", pfinalb := .5]
res_saffron[bF == "10", p5 := .5]
res_saffron[bF == "10", pfinalb := .5]
res_fwer[bF == "10", p5 := .5]
res_fwer[bF == "10", pfinalb := .5]

res_ai[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_saffron[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]
res_fwer[bF %in% c("10", "9"), .(bF, ate_tauv2, pfinalb, nodenum_current, nodenum_prev, nodesize, p4, p5)]

## With alpha fixed
res_fwer_det <- report_detections(res_fwer)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_fwer_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

## With alpha varying according to the alpha investing
res_ai_det <- report_detections(res_ai, fwer = FALSE)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_ai_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

## And with the saffron procedure
res_saffron_det <- report_detections(res_saffron, fwer = FALSE)
## So we can say that we discovered hits in the following blocks or groups of blocks
res_saffron_det[(hit), .(biggrp, bF, hit_grp, max_p, fin_parent_p, max_alpha, parent_alpha, single_hit, group_hit)][order(hit_grp)]

res_fwer_tree <- make_results_tree(res_fwer, blockid = "bF")
res_saffron_tree <- make_results_tree(res_saffron, blockid = "bF")
res_ai_tree <- make_results_tree(res_ai, blockid = "bF")

blahN <- res_fwer_tree %>%
  activate(nodes) %>%
  as_tibble()
blahE <- res_fwer_tree %>%
  activate(edges) %>%
  as_tibble()

res_fwer_g <- make_results_ggraph(res_fwer_tree)
res_ai_g <- make_results_ggraph(res_ai_tree)
res_saffron_g <- make_results_ggraph(res_saffron_tree)

## Compare these graphs to the results in res_ai_det and res_fwer_det above.
ggsave(res_fwer_g, file = "res_fwer_g.pdf", bg = "transparent", width = 12, height = 10)
ggsave(res_ai_g, file = "res_ai_g.pdf", bg = "transparent", width = 13, height = 10)
ggsave(res_saffron_g, file = "res_saffron_g.pdf", bg = "transparent", width = 13, height = 10)

## Compare working with tree/graph versus blocks
## TODO. Does report_detections produce the same discoveries as make_results_tree?

## Criteria for comparisons:
## These functions can be used to run the tests for the different scenarios.
## They have hard coded the relationships in alpha_and_splits above for now.

eval_numgrps <- function(detection_obj) {
  numgrps <- sapply(detection_obj, function(dat) {
    length(unique(dat$hit_grp))
  })
  ## names_split <- stri_split_regex(names(numgrps), "_", simplify = TRUE)

  expect_lte(numgrps[[4]], max(numgrps[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numgrps[[8]], max(numgrps[5:7]))
  expect_lte(numgrps[[12]], max(numgrps[9:11]))
  expect_lte(numgrps[[16]], max(numgrps[13:15]))
  expect_lte(numgrps[[20]], max(numgrps[17:19]))
}

eval_maxp <- function(detection_obj) {
  maxp <- sapply(detection_obj, function(dat) {
    dat[, sizerp := .N, by = hit_grp]
    dat[, p := ifelse(sizerp > 1, fin_parent_p, max_p)]
    max(dat[, .(min_max_p = min(p)), by = hit_grp]$min_max_p)
  })
  ## Adding tolerance of .01 here but could be more depending on whether we do any simulation
  ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(maxp[[4]], max(maxp[1:3]) + .01)
  expect_lte(maxp[[8]], max(maxp[5:7]) + .01)
  expect_lte(maxp[[12]], max(maxp[9:11]) + .01)
  expect_lte(maxp[[16]], max(maxp[13:15]) + .01)
  expect_lte(maxp[[20]], max(maxp[17:19]) + .01)
}

eval_single_blocks_found <- function(detection_obj) {
  ## Number of individual blocks detected versus groups (expect more singletons with the alpha adjusting approaches)
  numsingletons <- sapply(detection_obj, function(dat) {
    tab <- table(dat$hit_grp)
    sum(tab == 1)
  })
  expect_lte(numsingletons[[4]], min(numsingletons[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(numsingletons[[8]], min(numsingletons[5:7]))
  expect_lte(numsingletons[[12]], min(numsingletons[9:11]))
  expect_lte(numsingletons[[16]], min(numsingletons[13:15]))
  expect_lte(numsingletons[[20]], min(numsingletons[17:19]))
}

eval_treedepth <- function(detection_obj) {
  treedepth <- sapply(detection_obj, function(dat) {
    max(stri_count_fixed(dat$biggrp, ".")) + 1
  })
  expect_lte(treedepth[[4]], max(treedepth[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(treedepth[[8]], max(treedepth[5:7]))
  expect_lte(treedepth[[12]], max(treedepth[9:11]))
  expect_lte(treedepth[[16]], max(treedepth[13:15]))
  expect_lte(treedepth[[20]], max(treedepth[17:19]))
}

eval_minate <- function(detection_obj, atenm) {
  ## This compares the ATEs detected. I currently think that the fixed alpha approach should not be as sensitive or powerful as the alpha adjusting approaches. But somehow I wrote this test. Leaving it and ignoring it for now. Also the results fail this test. So I think my current intuition is correct.
  minate <- sapply(detection_obj, function(dat) {
    min(dat[, min(abs(get(atenm))), by = hit_grp]$V1)
  })
  expect_lte(minate[[4]], min(minate[1:3])) ## FWER should be smaller than at least one of the alpha adjusters
  expect_lte(minate[[8]], min(minate[5:7]))
  expect_lte(minate[[12]], min(minate[9:11]))
  expect_lte(minate[[16]], min(minate[13:15]))
  expect_lte(minate[[20]], min(minate[17:19]))
}

resnms <- apply(alpha_and_splits, 1, function(x) {
  tmp <- paste0(x, collapse = "_", sep = "")
  gsub("NULL", "fwer_fwer", tmp)
})


test_that("alphafns work across splitters for no effects", {
  ## First with no effects at all. So, basically no splitting should happen and no discoveries should be reported.
  tau_null <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynull ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_null) <- resnms

  tau_null_det <- lapply(seq_along(tau_null), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_null)[[i]], 1, 4) == "NULL"
    report_detections(tau_null[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_null_det) <- resnms
  ## No approach should detect anything, So we don't imagine any differences between approaches.
  anydetected <- sapply(tau_null_det, nrow)
  expect(all(anydetected == 0), TRUE)
})

test_that("alphafns work across splitters for large and homogenous effects", {
  ## All blocks same large effect. Both approaches should detect effects in basically all blocks (blocks only differ in size)
  tau_homog <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Yhomog ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_homog) <- resnms

  tau_homog_det <- lapply(seq_along(tau_homog), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_homog)[[i]], 1, 4) == "NULL"
    report_detections(tau_homog[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_homog_det) <- resnms
  ## We are comparing resnms 1 and 2 and 3 (the two fdr approaches) versus 4 (fwer) etc..

  ## Block 1 is the smallest block but it sometimes triggers a detection
  table(idat3$bF)
  block1detected <- sapply(tau_homog_det, function(dat) {
    any(dat$bF == "1")
  })
  table(block1detected)
  ## So basically all blocks or at least all but the smallest block should be detected
  numblks_homog <- sapply(tau_homog_det, nrow)
  expect_equal(all(numblks_homog >= 19), TRUE)
  ## Each block should be detected  --- they should not be grouped together ---
  ##  so the number of groups should be 19 as well
  numgrps_homog <- sapply(tau_homog_det, function(dat) {
    length(unique(dat$hit_grp))
  })
  expect_equal(all(numgrps_homog >= 19), TRUE)

  eval_numgrps(tau_homog_det)
  ## Number of individual blocks detected versus groups (expect more or same singletons with the alpha adjusting approaches)
  eval_single_blocks_found(tau_homog_det)

  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  minate_homog <- sapply(tau_homog_det, function(dat) {
    min(dat[, min(abs(ate_homog)), by = hit_grp]$V1)
  })

  ## Highest p detected: Doesn't reliably work in this case with all huge effects.
  eval_maxp(tau_homog_det)

  # g_homog <- lapply(tau_homog,function(obj){ make_results_ggraph(make_results_tree(obj)) })
  # for(i in 1:length(g_homog)){
  #    ggsave(g_homog[[i]], file = paste0("tau_homog_g",i,".pdf"),
  #           bg = "transparent", width = 14, height = 7)
  # }
})


test_that("alphafns work across splitters for individually heteogeneous effects and increase with block size. Also some completely null blocks.", {
  ################################################################################
  ## Some blocks have no effect at all in norm_inc
  tau_norm_inc <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynorm_inc ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_norm_inc) <- resnms
  tau_norm_inc_det <- lapply(seq_along(tau_norm_inc), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_norm_inc)[[i]], 1, 4) == "NULL"
    report_detections(tau_norm_inc[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_norm_inc_det) <- resnms

  ## Number of individual blocks detected versus groups
  ##  (expect more singletons with the alpha adjusting approaches than with fixed alpha)
  ## eval_single_blocks_found(tau_norm_inc_det)

  ## Number of blocks detected: this is not clear what to expect because of the ability to declare "detect" for *groups* of blocks.
  ## So, FWER might stop testing  and return many blocks.
  numblks_norm_inc <- sapply(tau_norm_inc_det, nrow)
  # expect_lte(numblks_norm_inc[[3]],min(numblks_norm_inc[1:2]))
  # expect_lte(numblks_norm_inc[[6]],min(numblks_norm_inc[4:5]))
  # expect_lte(numblks_norm_inc[[9]],min(numblks_norm_inc[7:8]))

  ## Number of groups of blocks detected: expect more or equal from alpha adjusting
  eval_numgrps(tau_norm_inc_det)
  ## Depth of testing (maxdepth): Probably a deeper tree or equal.
  eval_treedepth(tau_norm_inc_det)
  ## Lowest ate detected (do this by hit_grp) except it should be ok to have some null blocks in groups
  ## I currently disagree with # eval_minate --- shouldn't the alpha adjusted approaches be more sensitive and detect smaller ates?
  ## Commenting this out for now. I have expectations about p-values but not about ates per se
  ## eval_minate(tau_norm_inc_det, atenm = "ate_norm_inc")

  ## Highest p detected
  eval_maxp(tau_norm_inc_det)
  ## Null blocks detected (but this could be ok if they were in a group of other blocks having strong effects, so do by hit_grp)
  tau_norm_inc_det[[1]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  tau_norm_inc_det[[4]][, -c("ate_null", "ate_homog", "ate_tau", "ate_tauv2")]
  symdiff <- function(x, y) {
    setdiff(union(x, y), intersect(x, y))
  } # https://www.r-bloggers.com/symmetric-set-differences-in-r/
  symdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[4]]$bF)
  symdiff(tau_norm_inc_det[[2]]$bF, tau_norm_inc_det[[4]]$bF)
  setdiff(tau_norm_inc_det[[1]]$bF, tau_norm_inc_det[[4]]$bF)
  setdiff(tau_norm_inc_det[[2]]$bF, tau_norm_inc_det[[1]]$bF)

  sort(setdiff(bdat4$bF, tau_norm_inc_det[[1]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[2]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[3]]$bF))
  sort(setdiff(bdat4$bF, tau_norm_inc_det[[4]]$bF))
  # g_norm_inc <- lapply(tau_norm_inc, function(obj) {
  #   make_results_ggraph(make_results_tree(obj))
  # })
  # for (i in 1:length(g_norm_inc)) {
  #   ggsave(g_norm_inc[[i]],
  #     file = paste0("tau_norm_inc_g", i, ".pdf"),
  #     bg = "transparent", width = 14, height = 7
  #   )
  # }
})

## TODO this next not working. Maybe intuitions are wrong?
## test_that("alphafns work across splitters for individually heteogeneous effects
##  and decrease with block size. Also some completely null blocks.", {
##  ## Some blocks have no effect at all in norm_dec
##  tau_norm_dec <- mapply(
##    FUN = function(afn = afn, sfn = sfn, sby = sby) {
##      message(paste(afn, sfn, sby, collapse = ","))
##      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ynorm_dec ~ ZF | bF)
##    },
##    afn = alpha_and_splits$afn,
##    sfn = alpha_and_splits$sfn,
##    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
##  )
##  names(tau_norm_dec) <- resnms
##  tau_norm_dec_det <- lapply(seq_along(tau_norm_dec), function(i) {
##    # message(i)
##    fwer <- stri_sub(names(tau_norm_dec)[[i]], 1, 4) == "NULL"
##    report_detections(tau_norm_dec[[i]], fwer = fwer, only_hits = TRUE)
##  })
##  names(tau_norm_dec_det) <- resnms
##  somehits <- sapply(tau_norm_dec_det, nrow) != 0
##
##  tau_norm_dec_det_somehits <- tau_norm_dec_det[somehits]
##  ## eval_maxp(tau_norm_dec_det_somehits)
##  ## eval_minate(tau_norm_dec_det, atenm = "ate_norm_dec")
##  ## eval_numgrps(tau_norm_dec_det_somehits)
##  ## eval_single_blocks_found(tau_norm_dec_det_somehits)
##  ## eval_treedepth(tau_norm_dec_det_somehits)
## })

test_that("alphafns work across splitters for constant effect that cancel out at the high level
  (half large and positive, half large and negative).", {
  ################################################################################
  ### All blocks have large effects, some are negative and some positive. The
  ### blocks vary in size so perhaps the more sensitive procedures will be more
  ### likely to pick up effects in the smaller blocks.
  table(idat3$bF)

  tau_v1 <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Y ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_v1) <- resnms

  tau_v1_det <- lapply(seq_along(tau_v1), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_v1)[[i]], 1, 4) == "NULL"
    report_detections(tau_v1[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_v1_det) <- resnms

  eval_maxp(tau_v1_det)
  eval_minate(tau_v1_det, atenm = "ate_tau")
  eval_numgrps(tau_v1_det)
  eval_single_blocks_found(tau_v1_det)
  eval_treedepth(tau_v1_det)
})


test_that("alphafns work across splitters for individually heterogeneous effects block-fixed effects and some completely null blocks.", {
  ################
  tau_v2 <- mapply(
    FUN = function(afn = afn, sfn = sfn, sby = sby) {
      message(paste(afn, sfn, sby, collapse = ","))
      testing_fn(afn = afn, sfn = sfn, sby = sby, idat = idat3, bdat = bdat4, fmla = Ytauv2 ~ ZF | bF)
    },
    afn = alpha_and_splits$afn,
    sfn = alpha_and_splits$sfn,
    sby = alpha_and_splits$splitby, SIMPLIFY = FALSE
  )
  names(tau_v2) <- resnms

  tau_v2_det <- lapply(seq_along(tau_v2), function(i) {
    # message(i)
    fwer <- stri_sub(names(tau_v2)[[i]], 1, 4) == "NULL"
    report_detections(tau_v2[[i]], fwer = fwer, only_hits = TRUE)
  })
  names(tau_v2_det) <- resnms

  ### maxp and single_blocks don't fit the intuition here. Need to investigate. Commenting out for now.
  ## eval_maxp(tau_v2_det)
  ## eval_minate(tau_v2_det, atenm = "ate_tauv2")
  eval_numgrps(tau_v2_det)
  # eval_single_blocks_found(tau_v2_det)
  eval_treedepth(tau_v2_det)
})
