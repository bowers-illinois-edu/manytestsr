library(devtools)
library(data.table)
library(dplyr)
load_all()
source("tests/testthat/make_test_data.R")

# Run the same tests as in the failing case
set.seed(12345)
res_half <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL,
  local_adj_p_fn = NULL, simthresh = 20, sims = 1000, maxtest = 2000,
  thealpha = .05, thew0 = 0.05 - 0.001,
  fmla = Y_half_tau1 ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = FALSE, copydts = TRUE, stop_splitby_constant = TRUE,
  return_what = c("blocks", "nodes")
)

set.seed(12345)
res_half_hommel <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL,
  local_adj_p_fn = local_hommel_all_ps, simthresh = 20, sims = 1000, maxtest = 2000,
  thealpha = .05, thew0 = 0.05 - 0.001,
  fmla = Y_half_tau1 ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = FALSE, copydts = TRUE, stop_splitby_constant = TRUE,
  return_what = c("blocks", "nodes")
)

# Look at differences in block data structure
bdat_orig <- res_half$bdat
bdat_hommel <- res_half_hommel$bdat

cat("Original block data p-values:\n")
p_cols <- grep("^p[0-9]", names(bdat_orig), value = TRUE)
orig_p_summary <- bdat_orig[, lapply(.SD, function(x) length(unique(x[!is.na(x)]))), .SDcols = p_cols]
print(orig_p_summary)

cat("\nHommel block data p-values:\n") 
hommel_p_summary <- bdat_hommel[, lapply(.SD, function(x) length(unique(x[!is.na(x)]))), .SDcols = p_cols]
print(hommel_p_summary)

cat("\nComparing p2 values:\n")
cat("Original unique p2:", length(unique(bdat_orig$p2[!is.na(bdat_orig$p2)])), "\n")
cat("Hommel unique p2:", length(unique(bdat_hommel$p2[!is.na(bdat_hommel$p2)])), "\n")

cat("\nOriginal p2 values:\n")
print(unique(bdat_orig$p2[!is.na(bdat_orig$p2)]))
cat("\nHommel p2 values:\n")
print(unique(bdat_hommel$p2[!is.na(bdat_hommel$p2)]))

# Check if some p2 values are missing entirely in Hommel
orig_p2_unique <- sort(unique(bdat_orig$p2[!is.na(bdat_orig$p2)]))
hommel_p2_unique <- sort(unique(bdat_hommel$p2[!is.na(bdat_hommel$p2)]))

cat("\nMissing p2 values in Hommel:\n")
missing_p2 <- setdiff(orig_p2_unique, hommel_p2_unique)
print(missing_p2)

cat("\nExtra p2 values in Hommel:\n") 
extra_p2 <- setdiff(hommel_p2_unique, orig_p2_unique)
print(extra_p2)

# Check how many blocks have each p2 value
cat("\nOriginal p2 block counts:\n")
orig_p2_counts <- bdat_orig[!is.na(p2), .N, by = p2][order(p2)]
print(orig_p2_counts)

cat("\nHommel p2 block counts:\n")
hommel_p2_counts <- bdat_hommel[!is.na(p2), .N, by = p2][order(p2)]
print(hommel_p2_counts)

# Look for missing blocks in any p-value columns
for (pcol in p_cols) {
  orig_blocks_with_p <- nrow(bdat_orig[!is.na(get(pcol))])
  hommel_blocks_with_p <- nrow(bdat_hommel[!is.na(get(pcol))])
  if (orig_blocks_with_p != hommel_blocks_with_p) {
    cat("\nMismatch in", pcol, "- Original:", orig_blocks_with_p, "Hommel:", hommel_blocks_with_p, "\n")
  }
}