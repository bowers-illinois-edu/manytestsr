# data-raw/generate-toy-tree.R
# Re-generate the package dataset(s).
# Requires: TreeTestSim (dev-only)
# stopifnot(requireNamespace("TreeTestSim", quietly = TRUE))
library(TreeTestSim)
set.seed(20251108)
# If RNG stability matters across R versions:
# RNGversion("4.1.0")  # freeze algorithms if needed

##############################
#### Now generate complete k-ary tree with known null and non-null blocks
### (blocks with known no-effects or known to have non-zero effects)

## We think this kind of structure will cause trouble to unadjusted approaches.
## A tree with 4 nodes per levels and 3 levels after the root where 50% of the
## leaves are non-null the tree is defined at the level of the block We are not
## going to pass this vector below when we assess weak FWER control

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
testthat::expect_equal(parent_lookup[22][[1]], floor((22 - 2) / 4) + 1)
testthat::expect_equal(parent_lookup[22][[1]], 6)

# for a given node, return the vector of its ancestors+itself, *excluding* the
# NA‐parent root. "node" is an integer that counts from 1 at the root upwards
# in a structured order (see the math above about how to calculate the node
# numbers of parents, similar calculations are used for children in
# generate_tree_DT)

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

#  apply get_path to every node, paste the names with ".", and coerce to factor since
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
library(manytestsr)
bdt1[, split1 := splitSpecifiedFactorMulti(bid = node, x = lvls_fac)]
bdt1[, split2 := splitSpecifiedFactorMulti(node, lvls_fac), by = split1]
bdt1[, split3 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2)]

with(bdt1, table(split1, exclude = c()))
with(bdt1, ftable(split1 = split1, split2 = split2, split3 = split3))

library(dplyr)
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

## This next returns all 0s when splitting is no longer possible. So, consider
## watching out for this in the find_blocks code

bdt1[, split4 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3)]
bdt1[, split5 := splitSpecifiedFactorMulti(node, lvls_fac), by = interaction(split1, split2, split3, split4)]

testthat::expect_equal(as.numeric(as.character(unique(bdt1$split4))), 0)
testthat::expect_equal(as.numeric(as.character(unique(bdt1$split5))), 0)

bdt1[, c("split1", "split2", "split3", "split4", "split5") := NULL]

#### Done assessing the splitSpecifiedFactorMulti function on the data created by generate_tree_DT

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

## test using some non-null effect blocks so that the algorithm descends more into the tree
## This creates effects block-by-block.
idt[, y1_half_tau1 := create_effects(
  idat = idt, ybase = "y0", blockid = "bF", tau_fn = tau_norm, tau_size = 1, prop_blocks_0 = .5,
  non_null_blocks = "nonnull"
)]
## Check that the blocks with no effects actually have no effects
test_null_effects <- idt[nonnull == FALSE, sum(y1_half_tau1 - y0)]
testthat::expect_equal(unique(test_null_effects), 0)

test_not_null <- idt[nonnull == TRUE, sum(y1_half_tau1 - y0)]
testthat::expect_gt(abs(unique(test_not_null)), 0)
## Make an observed outcome
idt[, Y_half_tau1 := trt * y1_half_tau1 + (1 - trt) * y0]


# attr(toy_data, "data_version") <- "1.0.0"

usethis::use_data(idt, bdt1,
  overwrite = TRUE, compress = "bzip2", version = 3
)
