library(devtools)
library(data.table)
library(dplyr)
load_all()
source("tests/testthat/make_test_data.R")

# Run the exact same test as in test_adaptations.R
set.seed(12345)
res_half <- find_blocks(
  idat = idt, bdat = bdt1, blockid = "bF",
  splitfn = splitSpecifiedFactorMulti, pfn = pOneway, alphafn = NULL,
  local_adj_p_fn = NULL, simthresh = 20, sims = 1000, maxtest = 2000,
  thealpha = .05, thew0 = 0.05 - 0.001,
  fmla = Y_half_tau1 ~ trtF | bF, splitby = "lvls_fac",
  blocksize = "nb", trace = TRUE, copydts = TRUE, stop_splitby_constant = TRUE,
  return_what = c("blocks", "nodes")
)

# Original node_dat from find_blocks
orig_node_dat <- res_half$node_dat
cat("Original node_dat structure:\n")
print(orig_node_dat[order(depth, nodenum)])
cat("\nOriginal nodes by depth:\n")
print(orig_node_dat[, .N, by = depth])
cat("Total original nodes with non-NA p:", sum(!is.na(orig_node_dat$p)), "\n")

# Reconstructed from make_results_tree
res_half_tree <- make_results_tree(res_half$bdat, block_id = "bF")
recon_nodes <- res_half_tree$nodes
cat("\nReconstructed node structure:\n") 
print(recon_nodes[order(depth, name)])
cat("\nReconstructed nodes by depth:\n")
print(recon_nodes[, .N, by = depth])
cat("Total reconstructed nodes with non-NA p:", sum(!is.na(recon_nodes$p)), "\n")

# Compare the structures
cat("\nComparing original vs reconstructed:\n")
cat("Original has", nrow(orig_node_dat), "nodes, reconstructed has", nrow(recon_nodes), "nodes\n")

# Look at what nodes are missing
orig_nodes <- unique(orig_node_dat$nodenum)
recon_nodes_ids <- unique(recon_nodes$name)
missing_from_recon <- setdiff(orig_nodes, recon_nodes_ids)
extra_in_recon <- setdiff(recon_nodes_ids, orig_nodes)

cat("Missing from reconstruction:", length(missing_from_recon), "nodes:", paste(sort(missing_from_recon), collapse=", "), "\n")
cat("Extra in reconstruction:", length(extra_in_recon), "nodes:", paste(sort(extra_in_recon), collapse=", "), "\n")

# Look at block data to understand better
bdat <- res_half$bdat
cat("\nBlock data structure:\n")
print(bdat[, .(nodenum_current, nodenum_prev, biggrp, p1, p2, p3, p4)][1:20])

# Specifically examine depth 2 nodes
depth2_orig <- orig_node_dat[depth == 2]
cat("\nOriginal depth 2 nodes:\n")
print(depth2_orig[order(nodenum)])

# Look at p2 values in block data
p2_blocks <- bdat[!is.na(p2)]
cat("\nBlocks with p2 values:\n")
print(unique(p2_blocks[, .(nodenum_current, nodenum_prev, p2)])[order(p2)])

# Check what p2 values should map to which nodes
cat("\nUnique p2 values and their block mappings:\n")
p2_summary <- p2_blocks[, .(
  nodes_current = paste(sort(unique(nodenum_current)), collapse=","),
  nodes_prev = paste(sort(unique(nodenum_prev)), collapse=","),
  count = .N
), by = p2]
print(p2_summary[order(p2)])