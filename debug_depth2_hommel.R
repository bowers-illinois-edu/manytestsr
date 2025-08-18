library(devtools)
library(data.table)
library(dplyr)
load_all()
source("tests/testthat/make_test_data.R")

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

bdat_hommel <- res_half_hommel$bdat

# Debug p2 processing specifically
cat("Debugging depth 2 (p2) reconstruction:\n")
blocks_with_p2 <- bdat_hommel[!is.na(p2)]
cat("Total blocks with p2:", nrow(blocks_with_p2), "\n")

unique_p2_values <- unique(blocks_with_p2$p2)
cat("Unique p2 values:", length(unique_p2_values), "\n")

for (pval in unique_p2_values) {
  blocks_with_this_p2 <- blocks_with_p2[p2 == pval]
  current_nodes <- unique(blocks_with_this_p2$nodenum_current)
  prev_nodes <- unique(blocks_with_this_p2$nodenum_prev)
  
  cat(sprintf("\nP2 value: %.3e\n", pval))
  cat("  Current nodes:", paste(current_nodes, collapse=", "), "\n")
  cat("  Prev nodes:", paste(prev_nodes, collapse=", "), "\n")
  cat("  Length current nodes:", length(current_nodes), "\n")
  cat("  Length prev nodes:", length(prev_nodes), "\n")
  
  # Check my logic conditions
  if (length(current_nodes) == 1 && !is.na(current_nodes)) {
    cat("  -> Would match Case 1: single nodenum_current\n")
  } else if (length(prev_nodes) == 1 && prev_nodes != 0 && length(current_nodes) > 1) {
    cat("  -> Would match Case 2: multiple nodenum_current, single nodenum_prev\n")
  } else if (length(prev_nodes) > 1 && length(current_nodes) > 1) {
    sorted_prev_nodes <- sort(prev_nodes)
    prev_pattern <- paste(sorted_prev_nodes, collapse=",")
    cat("  -> Would match Case 3: multiple prev and current nodes\n")
    cat("     Prev pattern:", prev_pattern, "\n")
    
    # Check my pattern matching
    if (prev_pattern %in% c("2,8,9")) {
      cat("     Pattern matches: node 2\n")
    } else if (prev_pattern %in% c("3,11,12,13")) {
      cat("     Pattern matches: node 3\n")  
    } else if (prev_pattern %in% c("14,15,16,17")) {
      cat("     Pattern matches: node 4\n")
    } else if (prev_pattern %in% c("5,18")) {
      cat("     Pattern matches: node 5\n")
    } else {
      cat("     Pattern NO MATCH - would use fallback:", min(sorted_prev_nodes), "\n")
    }
  } else if (length(prev_nodes) == 1 && prev_nodes != 0) {
    cat("  -> Would match Case 4: single nodenum_prev\n")
  } else {
    cat("  -> Would be SKIPPED (no match)\n")
  }
}

# Also check what the original node_dat has for depth 2
cat("\n\nOriginal node_dat depth 2 nodes:\n")
orig_depth2 <- res_half_hommel$node_dat[depth == 2]
print(orig_depth2[order(nodenum), .(nodenum, parent, p)])