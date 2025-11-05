#' Intersection-Union Tests for Hierarchical Hypotheses
#'
#' Implementation of proper intersection-union testing for hierarchical structures.
#' This ensures that composite null hypotheses are tested appropriately and
#' that the hierarchical structure is respected in hypothesis testing.
#' Particularly useful for validating the logical structure of find_blocks() results.
#'
#' @param node_dat Node data with p-values and hierarchy
#' @param tracker Node tracker with tree structure
#' @param alpha Type I error rate
#' @param union_method Method for testing union alternatives ("max", "simes", "fisher")
#' @param intersection_method Method for testing intersection nulls ("min", "simes", "fisher")
#' @return Updated node data with intersection-union test results
#' @references
#' Berger, R. L. (1982). Multiparameter hypothesis testing and acceptance sampling.
#' Technometrics, 24(4), 295-300.
#' @examples
#' \dontrun{
#' # Apply intersection-union tests to find_blocks results
#' # Requires dplyr package
#' data(example_dat, package = "manytestsr")
#' library(data.table)
#' library(dplyr)
#' 
#' # Prepare data
#' idat <- as.data.table(example_dat)
#' bdat <- idat %>%
#'   group_by(blockF) %>%
#'   summarize(
#'     nb = n(),
#'     pb = mean(trt),
#'     hwt = (nb / nrow(idat)) * (pb * (1 - pb)),
#'     .groups = "drop"
#'   ) %>%
#'   as.data.table()
#' 
#' # Run find_blocks to get hierarchical structure
#' results <- find_blocks(
#'   idat = idat,
#'   bdat = bdat,
#'   blockid = "blockF",
#'   splitfn = splitCluster,
#'   pfn = pIndepDist,
#'   fmla = Y1 ~ trtF | blockF,
#'   splitby = "hwt",
#'   parallel = "no",
#'   maxtest = 15,
#'   trace = FALSE
#' )
#' 
#' # Apply intersection-union tests to validate hypothesis structure
#' iu_results <- intersection_union_tests(
#'   node_dat = results$node_dat,
#'   tracker = results$node_tracker,  # Assuming tracker is available
#'   alpha = 0.05,
#'   union_method = "simes",
#'   intersection_method = "simes"
#' )
#' 
#' # Examine intersection-union test results
#' if ("iu_p_intersection" %in% names(iu_results$node_dat)) {
#'   iu_summary <- iu_results$node_dat[, .(
#'     nodenum, 
#'     original_p = p,
#'     intersection_p = iu_p_intersection,
#'     union_p = iu_p_union,
#'     reject_intersection = iu_reject_intersection,
#'     reject_union = iu_reject_union
#'   )]
#'   
#'   cat("Intersection-Union Test Results:\n")
#'   print(head(iu_summary))
#'   
#'   # Count rejections by method
#'   cat("Intersection null rejections:", sum(iu_summary$reject_intersection, na.rm = TRUE), "\n")
#'   cat("Union alternative rejections:", sum(iu_summary$reject_union, na.rm = TRUE), "\n")
#' }
#' 
#' # Check consistency of intersection-union results
#' consistency_check <- check_iu_consistency(iu_results)
#' if (consistency_check$is_consistent) {
#'   cat("Intersection-union tests are logically consistent.\n")
#' } else {
#'   cat("Found inconsistencies:\n")
#'   print(consistency_check$inconsistencies)
#' }
#' 
#' # Visualize intersection-union results
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   iu_plot <- plot_intersection_union_results(iu_results)
#'   print(iu_plot)
#' }
#' }
#' @export
intersection_union_tests <- function(node_dat, tracker, alpha = 0.05,
                                   union_method = "max", 
                                   intersection_method = "simes") {
  
  # Identify the hierarchical structure
  hierarchy <- build_hierarchy_structure(node_dat, tracker)
  
  # For each node, define the intersection null and union alternative
  iu_results <- list()
  
  for (node_id in hierarchy$all_nodes) {
    # Get children of this node
    children <- hierarchy$children[[as.character(node_id)]]
    
    if (is.null(children) || length(children) == 0) {
      # Leaf node - simple hypothesis
      p_val <- node_dat[nodenum == node_id, p]
      iu_results[[as.character(node_id)]] <- list(
        node_id = node_id,
        is_leaf = TRUE,
        p_intersection = p_val,
        p_union = p_val,
        reject_intersection = !is.na(p_val) && p_val <= alpha,
        reject_union = !is.na(p_val) && p_val <= alpha
      )
    } else {
      # Internal node - composite hypothesis
      # H0: all children have null effects (intersection null)
      # H1: at least one child has non-null effect (union alternative)
      
      child_pvals <- sapply(children, function(cid) {
        p_val <- node_dat[nodenum == cid, p]
        if (length(p_val) == 0 || is.na(p_val)) return(1.0)
        return(p_val)
      })
      
      # Test intersection null H0: all children null
      p_intersection <- test_intersection_null(child_pvals, intersection_method)
      
      # Test union alternative H1: at least one child non-null
      p_union <- test_union_alternative(child_pvals, union_method)
      
      iu_results[[as.character(node_id)]] <- list(
        node_id = node_id,
        is_leaf = FALSE,
        children = children,
        child_pvals = child_pvals,
        p_intersection = p_intersection,
        p_union = p_union,
        reject_intersection = p_intersection <= alpha,
        reject_union = p_union <= alpha
      )
    }
  }
  
  # Add results to node_dat
  node_dat[, iu_p_intersection := NA_real_]
  node_dat[, iu_p_union := NA_real_]
  node_dat[, iu_reject_intersection := FALSE]
  node_dat[, iu_reject_union := FALSE]
  
  for (result in iu_results) {
    node_dat[nodenum == result$node_id, iu_p_intersection := result$p_intersection]
    node_dat[nodenum == result$node_id, iu_p_union := result$p_union] 
    node_dat[nodenum == result$node_id, iu_reject_intersection := result$reject_intersection]
    node_dat[nodenum == result$node_id, iu_reject_union := result$reject_union]
  }
  
  return(list(
    node_dat = node_dat,
    iu_results = iu_results,
    hierarchy = hierarchy
  ))
}

#' Build hierarchy structure from node data and tracker
#' @param node_dat Node data table
#' @param tracker Node tracker object
#' @return List with hierarchy information
build_hierarchy_structure <- function(node_dat, tracker) {
  
  all_nodes <- unique(node_dat$nodenum)
  children <- list()
  parents <- list()
  
  # Build parent-child relationships
  for (node_id in all_nodes) {
    node_children <- tracker$tracker[parent_id == node_id, node_id]
    node_parent <- tracker$tracker[node_id == node_id, parent_id]
    
    if (length(node_parent) > 0 && node_parent[1] != 0) {
      parents[[as.character(node_id)]] <- node_parent[1]
    }
    
    if (length(node_children) > 0) {
      children[[as.character(node_id)]] <- node_children
    }
  }
  
  # Identify levels
  max_depth <- max(node_dat$depth, na.rm = TRUE)
  levels <- list()
  for (d in 1:max_depth) {
    levels[[d]] <- node_dat[depth == d, nodenum]
  }
  
  return(list(
    all_nodes = all_nodes,
    children = children,
    parents = parents,
    levels = levels,
    max_depth = max_depth
  ))
}

#' Test intersection null hypothesis
#' @param p_values Vector of child p-values
#' @param method Method for combination
#' @return Combined p-value for intersection null
test_intersection_null <- function(p_values, method) {
  
  # Remove missing values
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0) return(1.0)
  
  switch(method,
    "min" = {
      # For intersection null, use minimum p-value (most liberal)
      # This corresponds to rejecting if ANY child hypothesis is rejected
      return(min(p_values))
    },
    "simes" = {
      # Simes method for intersection null
      k <- length(p_values)
      sorted_p <- sort(p_values)
      i_seq <- seq_len(k)
      simes_vals <- (k / i_seq) * sorted_p
      return(min(simes_vals))
    },
    "fisher" = {
      # Fisher's method
      if (any(p_values == 0)) return(0)
      chi_sq_stat <- -2 * sum(log(p_values))
      return(1 - pchisq(chi_sq_stat, df = 2 * length(p_values)))
    },
    {
      stop("Unknown intersection method: ", method)
    }
  )
}

#' Test union alternative hypothesis
#' @param p_values Vector of child p-values  
#' @param method Method for combination
#' @return Combined p-value for union alternative
test_union_alternative <- function(p_values, method) {
  
  # Remove missing values
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0) return(1.0)
  
  switch(method,
    "max" = {
      # For union alternative, use maximum p-value (most conservative) 
      # This corresponds to rejecting only if ALL child hypotheses are rejected
      return(max(p_values))
    },
    "simes" = {
      # Simes method (same as intersection but different interpretation)
      k <- length(p_values)
      sorted_p <- sort(p_values)
      i_seq <- seq_len(k)
      simes_vals <- (k / i_seq) * sorted_p
      return(min(simes_vals))
    },
    "fisher" = {
      # Fisher's method
      if (any(p_values == 0)) return(0)
      chi_sq_stat <- -2 * sum(log(p_values))
      return(1 - pchisq(chi_sq_stat, df = 2 * length(p_values)))
    },
    {
      stop("Unknown union method: ", method)
    }
  )
}

#' Check consistency of intersection-union test results
#' @param iu_results Results from intersection_union_tests
#' @return Logical vector indicating consistency
check_iu_consistency <- function(iu_results) {
  
  node_dat <- iu_results$node_dat
  hierarchy <- iu_results$hierarchy
  
  # Check logical consistency
  inconsistencies <- list()
  
  # Rule 1: If intersection null is rejected, union alternative should be rejected
  rule1_violations <- node_dat[
    iu_reject_intersection == TRUE & iu_reject_union == FALSE,
    nodenum
  ]
  
  if (length(rule1_violations) > 0) {
    inconsistencies$rule1 <- paste("Nodes", paste(rule1_violations, collapse = ", "), 
                                   "reject intersection null but not union alternative")
  }
  
  # Rule 2: For internal nodes, check parent-child consistency
  for (node_id in hierarchy$all_nodes) {
    children <- hierarchy$children[[as.character(node_id)]]
    if (!is.null(children) && length(children) > 0) {
      
      parent_rejects <- node_dat[nodenum == node_id, iu_reject_intersection]
      child_rejections <- node_dat[nodenum %in% children, iu_reject_intersection]
      
      # If parent rejects intersection null, at least one child should too
      if (length(parent_rejects) > 0 && parent_rejects && !any(child_rejections, na.rm = TRUE)) {
        inconsistencies[[paste0("node_", node_id)]] <- paste(
          "Node", node_id, "rejects intersection null but no children do"
        )
      }
    }
  }
  
  return(list(
    is_consistent = length(inconsistencies) == 0,
    inconsistencies = inconsistencies
  ))
}

#' Visualize intersection-union test results
#' @param iu_results Results from intersection_union_tests
#' @param highlight_inconsistencies Logical, whether to highlight inconsistencies
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_hline geom_vline labs scale_size_manual theme_minimal theme element_text
#' @export
plot_intersection_union_results <- function(iu_results, highlight_inconsistencies = TRUE) {

  node_dat <- iu_results$node_dat
  
  # Check for inconsistencies if requested
  consistency_check <- if (highlight_inconsistencies) {
    check_iu_consistency(iu_results)
  } else {
    list(is_consistent = TRUE, inconsistencies = list())
  }
  
  # Prepare data for plotting
  plot_data <- node_dat[, .(
    nodenum,
    depth,
    p_intersection = iu_p_intersection,
    p_union = iu_p_union,
    reject_intersection = iu_reject_intersection,
    reject_union = iu_reject_union
  )]
  
  # Add consistency flag
  plot_data[, consistent := TRUE]
  if (!consistency_check$is_consistent) {
    # Mark inconsistent nodes
    inconsistent_nodes <- as.numeric(gsub(".*node_([0-9]+).*", "\\1", 
                                        names(consistency_check$inconsistencies)))
    plot_data[nodenum %in% inconsistent_nodes, consistent := FALSE]
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = p_intersection, y = p_union)) +
    geom_point(aes(color = factor(depth), 
                   shape = ifelse(consistent, "Consistent", "Inconsistent"),
                   size = ifelse(reject_intersection | reject_union, "Rejected", "Not Rejected")),
               alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      x = "P-value for Intersection Null",
      y = "P-value for Union Alternative", 
      title = "Intersection-Union Test Results",
      subtitle = "Red dashed lines indicate alpha = 0.05",
      color = "Tree Depth",
      shape = "Consistency",
      size = "Decision"
    ) +
    scale_size_manual(values = c("Rejected" = 3, "Not Rejected" = 1.5)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  return(p)
}