#' Meinshausen's Hierarchical Testing with Sequential Rejection Principle
#'
#' Implementation of Meinshausen's (2008) hierarchical testing of variable importance
#' enhanced with the sequential rejection principle of Goeman and Solari (2010).
#' This approach provides a unified framework for hierarchical testing that maintains
#' FWER control while leveraging the hierarchical structure for increased power.
#'
#' @param node_dat Data.table containing node-level test results with columns:
#'   nodenum, p, depth, nodesize, and optionally parent, testable
#' @param node_tracker Node tracking object containing the hierarchical structure  
#' @param alpha Global Type I error rate (default: 0.05)
#' @param method Method for combining p-values ("simes", "fisher", "bonferroni")
#' @param use_sequential Logical, whether to use sequential rejection principle
#' @return Updated node_dat with Meinshausen hierarchical testing results
#' @references 
#' Meinshausen, N. (2008). Hierarchical testing of variable importance. 
#' Biometrika 95, 265-278.
#' 
#' Goeman, J. J., & Solari, A. (2010). The sequential rejection principle of 
#' familywise error control. Annals of Statistics 38, 3782-3810.
#' @export
meinshausen_hierarchical_test <- function(node_dat, node_tracker, alpha = 0.05, 
                                        method = "simes", use_sequential = TRUE) {
  
  # Parameter validation
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  if (!method %in% c("simes", "fisher", "bonferroni")) {
    stop("method must be one of 'simes', 'fisher', 'bonferroni'")
  }
  
  # Ensure we have a copy to avoid modifying original
  node_dat <- data.table::copy(node_dat)
  
  # Initialize results columns
  node_dat[, meinshausen_reject := FALSE]
  node_dat[, meinshausen_adjusted_p := p]
  node_dat[, meinshausen_level_tested := FALSE]
  
  # Step 1: Test global null hypothesis (root node)
  root_nodes <- node_dat[depth == 1 | (is.na(depth) & nodenum == 1)]
  
  if (nrow(root_nodes) == 0) {
    warning("No root node found in node_dat")
    return(node_dat)
  }
  
  # Global test using all leaf nodes
  if (use_sequential) {
    # Apply sequential rejection principle starting from root
    node_dat <- apply_sequential_rejection_meinshausen(node_dat, node_tracker, alpha, method)
  } else {
    # Apply traditional Meinshausen procedure
    node_dat <- apply_traditional_meinshausen(node_dat, node_tracker, alpha, method)
  }
  
  return(node_dat)
}

#' Apply Sequential Rejection Principle (Goeman-Solari enhancement)
#' @param node_dat Node data
#' @param node_tracker Node tracking structure
#' @param alpha Error rate
#' @param method P-value combination method
#' @return Updated node data with sequential rejections
apply_sequential_rejection_meinshausen <- function(node_dat, node_tracker, alpha, method) {
  
  # Sort nodes by depth (breadth-first traversal)
  data.table::setorder(node_dat, depth, nodenum)
  
  # Track which nodes can be tested (ancestors rejected)
  testable_nodes <- node_dat[depth == 1, nodenum]  # Start with root
  
  # Sequential testing by level
  max_depth <- max(node_dat$depth, na.rm = TRUE)
  
  for (current_depth in 1:max_depth) {
    level_nodes <- node_dat[depth == current_depth & nodenum %in% testable_nodes]
    
    if (nrow(level_nodes) == 0) next
    
    # Apply level-specific adjustment
    level_adjusted_p <- adjust_pvalues_level(level_nodes$p, method, alpha)
    
    # Update adjusted p-values
    for (i in seq_len(nrow(level_nodes))) {
      node_idx <- which(node_dat$nodenum == level_nodes$nodenum[i])
      node_dat[node_idx, meinshausen_adjusted_p := level_adjusted_p[i]]
      node_dat[node_idx, meinshausen_level_tested := TRUE]
      
      # Apply rejection rule
      if (!is.na(level_adjusted_p[i]) && level_adjusted_p[i] <= alpha) {
        node_dat[node_idx, meinshausen_reject := TRUE]
        
        # Add children to testable nodes for next level
        children <- get_children_nodes(level_nodes$nodenum[i], node_tracker)
        testable_nodes <- c(testable_nodes, children)
      }
    }
  }
  
  return(node_dat)
}

#' Apply Traditional Meinshausen Procedure
#' @param node_dat Node data
#' @param node_tracker Node tracking structure  
#' @param alpha Error rate
#' @param method P-value combination method
#' @return Updated node data with traditional Meinshausen results
apply_traditional_meinshausen <- function(node_dat, node_tracker, alpha, method) {
  
  # Start with root node test
  root_nodes <- node_dat[depth == 1]
  
  if (nrow(root_nodes) == 0) return(node_dat)
  
  # Global test combining all p-values
  global_p <- combine_pvalues_traditional(node_dat$p, method)
  
  # Mark root as tested
  node_dat[depth == 1, meinshausen_level_tested := TRUE]
  node_dat[depth == 1, meinshausen_adjusted_p := global_p]
  
  # If global test rejects, proceed down hierarchy
  if (!is.na(global_p) && global_p <= alpha) {
    node_dat[depth == 1, meinshausen_reject := TRUE]
    
    # Recursively test children
    for (root_id in root_nodes$nodenum) {
      node_dat <- test_children_recursive_meinshausen(
        root_id, node_dat, node_tracker, alpha, method
      )
    }
  }
  
  return(node_dat)
}

#' Recursive testing of children in Meinshausen procedure
#' @param parent_id Parent node ID
#' @param node_dat Node data
#' @param node_tracker Node tracking structure
#' @param alpha Error rate  
#' @param method P-value combination method
#' @return Updated node data
test_children_recursive_meinshausen <- function(parent_id, node_dat, node_tracker, alpha, method) {
  
  # Get children of current parent
  children <- get_children_nodes(parent_id, node_tracker)
  
  if (length(children) == 0) return(node_dat)
  
  # Test each child cluster
  for (child_id in children) {
    child_nodes <- get_descendant_nodes(child_id, node_tracker, node_dat)
    
    if (nrow(child_nodes) == 0) next
    
    # Combine p-values for this cluster
    cluster_p <- combine_pvalues_traditional(child_nodes$p, method)
    
    # Adjust for multiple testing at this level
    n_siblings <- length(children)
    adjusted_p <- min(1, cluster_p * n_siblings)  # Bonferroni correction
    
    # Update node data
    child_idx <- which(node_dat$nodenum == child_id)
    if (length(child_idx) > 0) {
      node_dat[child_idx, meinshausen_level_tested := TRUE]
      node_dat[child_idx, meinshausen_adjusted_p := adjusted_p]
      
      if (!is.na(adjusted_p) && adjusted_p <= alpha) {
        node_dat[child_idx, meinshausen_reject := TRUE]
        
        # Recursively test children of this child
        node_dat <- test_children_recursive_meinshausen(
          child_id, node_dat, node_tracker, alpha, method
        )
      }
    }
  }
  
  return(node_dat)
}

#' Adjust p-values at a given level using sequential rejection
#' @param p_values Vector of p-values at current level
#' @param method Adjustment method
#' @param alpha Error rate
#' @return Vector of adjusted p-values
adjust_pvalues_level <- function(p_values, method, alpha) {
  
  n <- length(p_values)
  if (n == 0) return(numeric(0))
  if (n == 1) return(p_values)
  
  # Sort p-values for sequential testing
  sorted_idx <- order(p_values)
  sorted_p <- p_values[sorted_idx]
  
  # Apply method-specific adjustment
  adjusted_p <- switch(method,
    "simes" = {
      # Simes-type adjustment for sequential rejection
      sapply(seq_along(sorted_p), function(i) {
        min(1, sorted_p[i] * n / i)
      })
    },
    "bonferroni" = {
      # Conservative Bonferroni adjustment
      pmin(1, sorted_p * n)
    },
    "fisher" = {
      # Fisher combination with sequential adjustment
      # This is more complex and approximated here
      fisher_stats <- -2 * log(sorted_p)
      chi_sq_p <- pchisq(fisher_stats, df = 2, lower.tail = FALSE)
      pmin(1, chi_sq_p * (n:1))
    }
  )
  
  # Enforce monotonicity (required for sequential rejection)
  adjusted_p <- cummax(adjusted_p)
  
  # Restore original order
  result <- numeric(n)
  result[sorted_idx] <- adjusted_p
  
  return(result)
}

#' Get children nodes from tracker
#' @param parent_id Parent node ID
#' @param node_tracker Node tracking structure
#' @return Vector of child node IDs
get_children_nodes <- function(parent_id, node_tracker) {
  if (is.null(node_tracker$tracker)) return(integer(0))
  
  tracker_dt <- node_tracker$tracker
  child_mask <- tracker_dt$parent_id == parent_id
  children <- tracker_dt$node_id[child_mask]
  return(children[!is.na(children)])
}

#' Get all descendant nodes for a given node
#' @param node_id Node ID
#' @param node_tracker Node tracking structure
#' @param node_dat Node data
#' @return Data table of descendant nodes
get_descendant_nodes <- function(node_id, node_tracker, node_dat) {
  
  # Get all nodes that are descendants of node_id
  descendants <- c()
  to_check <- node_id
  max_iterations <- 100  # Safety limit to prevent infinite loops
  iteration <- 0
  
  while (length(to_check) > 0 && iteration < max_iterations) {
    current <- to_check[1]
    to_check <- to_check[-1]
    
    # Avoid duplicates
    if (!current %in% descendants) {
      descendants <- c(descendants, current)
      
      # Add children to check
      children <- get_children_nodes(current, node_tracker)
      # Only add children that aren't already in descendants to avoid cycles
      new_children <- children[!children %in% descendants]
      to_check <- c(to_check, new_children)
    }
    
    iteration <- iteration + 1
  }
  
  if (iteration >= max_iterations) {
    warning("Maximum iterations reached in get_descendant_nodes - possible infinite loop prevented")
  }
  
  # Return node data for descendants
  return(node_dat[nodenum %in% descendants])
}

#' Traditional p-value combination (used in original Meinshausen)
#' @param p_values Vector of p-values to combine
#' @param method Combination method
#' @return Combined p-value
combine_pvalues_traditional <- function(p_values, method) {
  
  if (length(p_values) == 0) return(1)
  
  # Remove NA values first  
  p_values <- p_values[!is.na(p_values)]
  if (length(p_values) == 0) return(1)
  if (length(p_values) == 1) return(p_values[1])
  
  switch(method,
    "simes" = {
      n <- length(p_values)
      sorted_p <- sort(p_values)
      min(sapply(seq_along(sorted_p), function(i) {
        min(1, sorted_p[i] * n / i)
      }))
    },
    "fisher" = {
      # Fisher's method
      if (any(p_values == 0)) return(0)
      fisher_stat <- -2 * sum(log(p_values))
      pchisq(fisher_stat, df = 2 * length(p_values), lower.tail = FALSE)
    },
    "bonferroni" = {
      min(1, min(p_values) * length(p_values))
    }
  )
}

#' Generate Meinshausen hierarchical clustering for variables
#'
#' Creates a hierarchical clustering structure suitable for Meinshausen testing
#' when applied to variable selection problems. This is typically used when
#' the goal is to test groups of correlated variables hierarchically.
#'
#' @param correlation_matrix Correlation matrix of variables
#' @param method Clustering method ("complete", "average", "single")
#' @param min_cluster_size Minimum size for clusters to be tested
#' @return List with clustering structure compatible with node_tracker format
#' @export
generate_meinshausen_hierarchy <- function(correlation_matrix, method = "complete", 
                                         min_cluster_size = 2) {
  
  if (!is.matrix(correlation_matrix)) {
    stop("correlation_matrix must be a matrix")
  }
  
  # Create distance matrix from correlation
  dist_matrix <- as.dist(1 - abs(correlation_matrix))
  
  # Hierarchical clustering
  hc <- hclust(dist_matrix, method = method)
  
  # Convert to node_tracker format
  n_vars <- nrow(correlation_matrix)
  n_nodes <- n_vars + nrow(hc$merge)
  
  tracker_data <- data.table::data.table(
    node_id = integer(),
    parent_id = integer(), 
    depth = integer(),
    cluster_size = integer(),
    variables = list()
  )
  
  # Add leaf nodes (individual variables)
  for (i in seq_len(n_vars)) {
    tracker_data <- rbind(tracker_data, data.table::data.table(
      node_id = i,
      parent_id = NA_integer_,
      depth = max(hc$height) + 1,
      cluster_size = 1L,
      variables = list(i)
    ))
  }
  
  # Add internal nodes from clustering
  for (i in seq_len(nrow(hc$merge))) {
    node_id <- n_vars + i
    left_child <- hc$merge[i, 1]
    right_child <- hc$merge[i, 2]
    
    # Convert negative indices to positive (leaf nodes)
    if (left_child < 0) left_child <- abs(left_child)
    else left_child <- left_child + n_vars
    
    if (right_child < 0) right_child <- abs(right_child)
    else right_child <- right_child + n_vars
    
    # Get variables in this cluster
    left_vars <- tracker_data[node_id == left_child, variables][[1]]
    right_vars <- tracker_data[node_id == right_child, variables][[1]]
    cluster_vars <- c(left_vars, right_vars)
    
    tracker_data <- rbind(tracker_data, data.table::data.table(
      node_id = node_id,
      parent_id = NA_integer_,
      depth = hc$height[i],
      cluster_size = length(cluster_vars),
      variables = list(cluster_vars)
    ))
    
    # Update parent information
    tracker_data[node_id == left_child, parent_id := node_id]
    tracker_data[node_id == right_child, parent_id := node_id]
  }
  
  # Convert heights to depth levels (reverse order)
  max_height <- max(tracker_data$depth, na.rm = TRUE)
  tracker_data[, depth := max_height - depth + 1]
  
  # Root node has no parent
  root_node <- tracker_data[depth == 1, node_id][1]
  tracker_data[node_id == root_node, parent_id := 0L]
  
  return(list(
    tracker = tracker_data,
    next_id = max(tracker_data$node_id) + 1L,
    hclust_object = hc,
    variable_names = rownames(correlation_matrix)
  ))
}

#' Integration function for find_blocks with Meinshausen testing
#'
#' Integrates Meinshausen hierarchical testing into the find_blocks workflow
#' 
#' @param node_dat Node-level data from find_blocks
#' @param node_tracker Node tracking structure  
#' @param alpha Error rate
#' @param method P-value combination method
#' @param use_sequential Use sequential rejection principle
#' @return Enhanced node_dat with Meinshausen results
#' @export
integrate_meinshausen_find_blocks <- function(node_dat, node_tracker, alpha = 0.05,
                                            method = "simes", use_sequential = TRUE) {
  
  tryCatch({
    # Apply Meinshausen hierarchical testing
    enhanced_node_dat <- meinshausen_hierarchical_test(
      node_dat = node_dat,
      node_tracker = node_tracker, 
      alpha = alpha,
      method = method,
      use_sequential = use_sequential
    )
    
    return(enhanced_node_dat)
    
  }, error = function(e) {
    warning("Meinshausen hierarchical testing failed: ", e$message, 
            ". Continuing with standard results.")
    return(node_dat)
  })
}