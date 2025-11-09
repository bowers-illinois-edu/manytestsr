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

# Compare the node structures
cat("Original (no adjustment) nodes:\n")
orig_nodes <- res_half$node_dat
print(orig_nodes[order(depth, nodenum), .(nodenum, parent, depth, p)])

cat("\nHommel adjusted nodes:\n")
hommel_nodes <- res_half_hommel$node_dat
print(hommel_nodes[order(depth, nodenum), .(nodenum, parent, depth, p)])

cat("\nNode count comparison:\n")
cat("Original:", nrow(orig_nodes), "nodes\n")
cat("Hommel:", nrow(hommel_nodes), "nodes\n")

# Get the reconstructed trees
res_half_tree <- make_results_tree(res_half, block_id = "bF")
res_half_hommel_tree <- make_results_tree(res_half_hommel, block_id = "bF")

cat("\nReconstructed node count comparison:\n")
cat("Original reconstructed:", nrow(res_half_tree$nodes), "nodes\n")
cat("Hommel reconstructed:", nrow(res_half_hommel_tree$nodes), "nodes\n")

cat("\nReconstructed nodes - Original:\n")
print(res_half_tree$nodes[order(depth, name), .(name, parent_name, depth, p)])

cat("\nReconstructed nodes - Hommel:\n")
print(res_half_hommel_tree$nodes[order(depth, name), .(name, parent_name, depth, p)])

# Compare specific columns that are failing in the tests
cat("\nComparing node_number columns:\n")
orig_node_numbers <- res_half_tree$nodes$node_number
hommel_node_numbers <- res_half_hommel_tree$nodes$node_number

cat("Original node_number length:", length(orig_node_numbers), "\n")
cat("Hommel node_number length:", length(hommel_node_numbers), "\n")

if (length(orig_node_numbers) != length(hommel_node_numbers)) {
  cat("Different lengths detected!\n")
  cat("Original node_numbers:", paste(orig_node_numbers, collapse = ", "), "\n")
  cat("Hommel node_numbers:", paste(hommel_node_numbers, collapse = ", "), "\n")
}
