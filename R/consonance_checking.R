#' Consonance Property Checking for Hierarchical Testing
#'
#' Implementation of consonance property validation for hierarchical multiple testing procedures.
#' The consonance property ensures logical consistency: if a hypothesis is rejected,
#' all hypotheses that logically imply it should also be rejectable.
#' Essential for validating find_blocks() results and ensuring logical coherence.
#'
#' @param node_dat Node data with test results
#' @param tracker Node tracker with hierarchy structure
#' @param rejection_column Column name containing rejection decisions
#' @param alpha Type I error rate used for testing
#' @return List with consonance check results and recommendations
#' @references
#' Goeman, J. J., & Solari, A. (2011). Multiple testing for exploratory research.
#' Statistical science, 26(4), 584-597.
#' @examples
#' \donttest{
#' # Check consonance of find_blocks results
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
#' # Run find_blocks
#' results <- find_blocks(
#'   idat = idat,
#'   bdat = bdat,
#'   blockid = "blockF",
#'   splitfn = splitCluster,
#'   pfn = pOneway,
#'   fmla = Y1 ~ trtF | blockF,
#'   splitby = "hwt",
#'   parallel = "no",
#'   maxtest = 20,
#'   trace = FALSE
#' )
#' 
#' # Check consonance of the testable decisions
#' consonance_check <- check_consonance_property(
#'   node_dat = results$node_dat,
#'   tracker = results$node_tracker,  # Assuming tracker is available
#'   rejection_column = "testable",
#'   alpha = 0.05
#' )
#' 
#' # Examine consonance results
#' if (consonance_check$is_consonant) {
#'   cat("The hierarchical testing procedure is consonant.\n")
#' } else {
#'   cat("Found consonance violations:\n")
#'   for (i in seq_along(consonance_check$violations)) {
#'     violation <- consonance_check$violations[[i]]
#'     cat("-", violation$violation_description, "\n")
#'   }
#'   
#'   # Print violation analysis
#'   analysis <- consonance_check$violation_analysis
#'   cat("\nViolation Analysis:\n")
#'   cat("Total violations:", analysis$total_violations, "\n")
#'   cat("Severity level:", analysis$severity, "\n")
#'   cat("Affected nodes:", paste(analysis$affected_nodes, collapse = ", "), "\n")
#' }
#' 
#' # Apply consonance correction if needed
#' if (!consonance_check$is_consonant) {
#'   cat("\nApplying consonance correction...\n")
#'   
#'   corrected_results <- apply_consonance_correction(
#'     node_dat = results$node_dat,
#'     tracker = results$node_tracker,
#'     rejection_column = "testable",
#'     correction_method = "propagate_up"
#'   )
#'   
#'   # Check if correction worked
#'   corrected_check <- check_consonance_property(
#'     corrected_results,
#'     results$node_tracker,
#'     "testable_consonant"
#'   )
#'   
#'   if (corrected_check$is_consonant) {
#'     cat("Consonance correction successful.\n")
#'   } else {
#'     cat("Additional corrections may be needed.\n")
#'   }
#' }
#' 
#' # Compare original and corrected rejection counts
#' original_rejections <- sum(results$node_dat$testable, na.rm = TRUE)
#' if (exists("corrected_results")) {
#'   corrected_rejections <- sum(corrected_results$testable_consonant, na.rm = TRUE)
#'   cat("Original rejections:", original_rejections, "\n")
#'   cat("Corrected rejections:", corrected_rejections, "\n")
#' }
#' }
#' @export
check_consonance_property <- function(node_dat, tracker, 
                                    rejection_column = "testable", 
                                    alpha = 0.05) {
  
  # Build hierarchical relationships
  hierarchy <- extract_hierarchy_relationships(node_dat, tracker)
  
  # Check consonance violations
  consonance_violations <- find_consonance_violations(
    node_dat, hierarchy, rejection_column
  )
  
  # Analyze violation patterns
  violation_analysis <- analyze_violation_patterns(consonance_violations, hierarchy)
  
  # Generate corrective actions
  corrective_actions <- generate_corrective_actions(
    consonance_violations, node_dat, hierarchy
  )
  
  # Validate logical hierarchy
  logical_consistency <- validate_logical_hierarchy(hierarchy, node_dat)
  
  return(list(
    is_consonant = length(consonance_violations) == 0,
    violations = consonance_violations,
    violation_analysis = violation_analysis,
    corrective_actions = corrective_actions,
    logical_consistency = logical_consistency,
    hierarchy = hierarchy
  ))
}

#' Extract hierarchy relationships
#' @param node_dat Node data table
#' @param tracker Node tracker object
#' @return List with hierarchy relationships
extract_hierarchy_relationships <- function(node_dat, tracker) {
  
  # Get all nodes and their relationships
  all_nodes <- unique(node_dat$nodenum)
  
  # Build ancestor and descendant relationships
  ancestors <- list()
  descendants <- list()
  direct_children <- list()
  direct_parent <- list()
  
  for (node_id in all_nodes) {
    # Find ancestors (path to root)
    ancestor_path <- build_numeric_ancestry(tracker, node_id)
    ancestors[[as.character(node_id)]] <- setdiff(ancestor_path, node_id)
    
    # Find descendants (all nodes below this one)
    node_descendants <- find_all_descendants(node_id, tracker, all_nodes)
    descendants[[as.character(node_id)]] <- node_descendants
    
    # Direct relationships
    children <- tracker$tracker[parent_id == node_id, node_id]
    if (length(children) > 0) {
      direct_children[[as.character(node_id)]] <- children
    }
    
    parent <- tracker$tracker[node_id == node_id, parent_id]
    if (length(parent) > 0 && parent[1] != 0) {
      direct_parent[[as.character(node_id)]] <- parent[1]
    }
  }
  
  return(list(
    all_nodes = all_nodes,
    ancestors = ancestors,
    descendants = descendants,
    direct_children = direct_children,
    direct_parent = direct_parent
  ))
}

#' Find all descendants of a node
#' @param node_id Node to find descendants for
#' @param tracker Node tracker object
#' @param all_nodes Vector of all node IDs
#' @return Vector of descendant node IDs
find_all_descendants <- function(node_id, tracker, all_nodes) {
  descendants <- c()
  
  # Use breadth-first search
  queue <- node_id
  visited <- c()
  
  while (length(queue) > 0) {
    current <- queue[1]
    queue <- queue[-1]
    
    if (current %in% visited) next
    visited <- c(visited, current)
    
    # Find direct children
    children <- tracker$tracker[parent_id == current, node_id]
    for (child in children) {
      if (!child %in% visited) {
        descendants <- c(descendants, child)
        queue <- c(queue, child)
      }
    }
  }
  
  return(descendants)
}

#' Find consonance violations
#' @param node_dat Node data table
#' @param hierarchy Hierarchy relationships
#' @param rejection_column Column with rejection decisions
#' @return List of violations
find_consonance_violations <- function(node_dat, hierarchy, rejection_column) {
  
  violations <- list()
  
  # Check each rejected hypothesis
  rejected_nodes <- node_dat[get(rejection_column) == TRUE, nodenum]
  
  for (rejected_node in rejected_nodes) {
    # Get ancestors of this rejected node
    node_ancestors <- hierarchy$ancestors[[as.character(rejected_node)]]
    
    if (length(node_ancestors) > 0) {
      # Check if all ancestors are also rejectable
      for (ancestor in node_ancestors) {
        ancestor_rejectable <- node_dat[nodenum == ancestor, get(rejection_column)]
        
        if (length(ancestor_rejectable) > 0 && !ancestor_rejectable) {
          # Consonance violation: descendant rejected but ancestor not rejectable
          violations[[length(violations) + 1]] <- list(
            type = "ancestor_not_rejectable",
            rejected_node = rejected_node,
            non_rejectable_ancestor = ancestor,
            violation_description = paste0(
              "Node ", rejected_node, " is rejected but ancestor ", 
              ancestor, " is not rejectable"
            )
          )
        }
      }
    }
  }
  
  # Check for upward consonance: if all children rejected, parent should be rejectable
  for (node_id in hierarchy$all_nodes) {
    children <- hierarchy$direct_children[[as.character(node_id)]]
    
    if (length(children) > 0) {
      child_rejections <- node_dat[nodenum %in% children, get(rejection_column)]
      
      # If all children are rejected
      if (length(child_rejections) > 0 && all(child_rejections, na.rm = TRUE)) {
        parent_rejectable <- node_dat[nodenum == node_id, get(rejection_column)]
        
        if (length(parent_rejectable) > 0 && !parent_rejectable) {
          violations[[length(violations) + 1]] <- list(
            type = "parent_not_rejectable_despite_all_children",
            non_rejectable_parent = node_id,
            rejected_children = children,
            violation_description = paste0(
              "All children (", paste(children, collapse = ", "), 
              ") of node ", node_id, " are rejected but parent is not rejectable"
            )
          )
        }
      }
    }
  }
  
  return(violations)
}

#' Analyze violation patterns
#' @param violations List of consonance violations
#' @param hierarchy Hierarchy structure
#' @return Analysis of violation patterns
analyze_violation_patterns <- function(violations, hierarchy) {
  
  if (length(violations) == 0) {
    return(list(
      total_violations = 0,
      violation_types = character(0),
      affected_nodes = numeric(0),
      severity = "none"
    ))
  }
  
  # Categorize violations
  violation_types <- sapply(violations, function(v) v$type)
  violation_counts <- table(violation_types)
  
  # Identify most problematic nodes
  affected_nodes <- unique(c(
    sapply(violations, function(v) v$rejected_node),
    sapply(violations, function(v) v$non_rejectable_ancestor),
    sapply(violations, function(v) v$non_rejectable_parent)
  ))
  affected_nodes <- affected_nodes[!is.na(affected_nodes)]
  
  # Assess severity
  severity <- if (length(violations) > length(hierarchy$all_nodes) * 0.1) {
    "high"
  } else if (length(violations) > 3) {
    "medium"
  } else {
    "low"
  }
  
  return(list(
    total_violations = length(violations),
    violation_types = violation_counts,
    affected_nodes = affected_nodes,
    severity = severity,
    violation_rate = length(violations) / length(hierarchy$all_nodes)
  ))
}

#' Generate corrective actions for consonance violations
#' @param violations List of violations
#' @param node_dat Node data table
#' @param hierarchy Hierarchy structure
#' @return Recommended corrective actions
generate_corrective_actions <- function(violations, node_dat, hierarchy) {
  
  if (length(violations) == 0) {
    return(list(message = "No corrective actions needed - consonance property satisfied"))
  }
  
  actions <- list()
  
  # Group violations by type
  violation_types <- sapply(violations, function(v) v$type)
  
  # Actions for ancestor not rejectable violations
  ancestor_violations <- violations[violation_types == "ancestor_not_rejectable"]
  if (length(ancestor_violations) > 0) {
    actions$ancestor_fixes <- list(
      action_type = "make_ancestors_rejectable",
      description = "Make ancestor nodes rejectable when descendants are rejected",
      affected_nodes = unique(sapply(ancestor_violations, function(v) v$non_rejectable_ancestor)),
      implementation = "Consider using proper closed testing procedure or adjust alpha levels upward in the hierarchy"
    )
  }
  
  # Actions for parent not rejectable violations
  parent_violations <- violations[violation_types == "parent_not_rejectable_despite_all_children"]
  if (length(parent_violations) > 0) {
    actions$parent_fixes <- list(
      action_type = "make_parents_rejectable",
      description = "Make parent nodes rejectable when all children are rejected",
      affected_nodes = unique(sapply(parent_violations, function(v) v$non_rejectable_parent)),
      implementation = "Apply intersection-union test principle or use appropriate combining function"
    )
  }
  
  # General recommendations
  actions$general_recommendations <- list(
    "Use proper closed testing procedure to ensure consonance",
    "Consider Bonferroni-Holm or other consonant procedures",
    "Validate that the hierarchical structure reflects logical relationships",
    "Check that p-value combination methods are appropriate for the hypothesis structure"
  )
  
  return(actions)
}

#' Validate logical hierarchy structure
#' @param hierarchy Hierarchy object
#' @param node_dat Node data table
#' @return Logical consistency check results
validate_logical_hierarchy <- function(hierarchy, node_dat) {
  
  issues <- list()
  
  # Check for cycles
  for (node_id in hierarchy$all_nodes) {
    ancestors <- hierarchy$ancestors[[as.character(node_id)]]
    if (node_id %in% ancestors) {
      issues$cycles <- c(issues$cycles, node_id)
    }
  }
  
  # Check for orphaned nodes (except root)
  # Create a logical vector for nodes without parents
  has_no_parent <- sapply(hierarchy$all_nodes, function(nid) {
    parent_info <- hierarchy$direct_parent[[as.character(nid)]]
    is.null(parent_info) || length(parent_info) == 0
  })
  
  roots <- hierarchy$all_nodes[has_no_parent]
  
  if (length(roots) > 1) {
    issues$multiple_roots <- roots
  }
  
  # Check depth consistency
  for (node_id in hierarchy$all_nodes) {
    node_depth <- node_dat[nodenum == node_id, depth]
    ancestors <- hierarchy$ancestors[[as.character(node_id)]]
    
    for (ancestor in ancestors) {
      ancestor_depth <- node_dat[nodenum == ancestor, depth]
      if (length(ancestor_depth) > 0 && length(node_depth) > 0 && 
          ancestor_depth >= node_depth) {
        issues$depth_inconsistency <- c(issues$depth_inconsistency, 
                                       paste0(ancestor, "->", node_id))
      }
    }
  }
  
  return(list(
    is_valid = length(issues) == 0,
    issues = issues,
    recommendations = if (length(issues) > 0) {
      c("Fix hierarchy structure before applying consonance checks",
        "Ensure proper parent-child relationships",
        "Validate tree construction algorithm")
    } else {
      "Hierarchy structure is logically consistent"
    }
  ))
}

#' Apply consonance correction to test results
#' @param node_dat Node data table with test results
#' @param tracker Node tracker object
#' @param rejection_column Column with rejection decisions
#' @param correction_method Method for correction ("propagate_up", "propagate_down", "closed_testing")
#' @return Corrected node data table
#' @export
apply_consonance_correction <- function(node_dat, tracker, 
                                      rejection_column = "testable",
                                      correction_method = "propagate_up") {
  
  hierarchy <- extract_hierarchy_relationships(node_dat, tracker)
  corrected_data <- copy(node_dat)
  
  # Create corrected rejection column
  corrected_column <- paste0(rejection_column, "_consonant")
  corrected_data[, (corrected_column) := get(rejection_column)]
  
  switch(correction_method,
    "propagate_up" = {
      # If a node is rejected, make all ancestors rejectable
      rejected_nodes <- corrected_data[get(rejection_column) == TRUE, nodenum]
      
      for (rejected_node in rejected_nodes) {
        ancestors <- hierarchy$ancestors[[as.character(rejected_node)]]
        for (ancestor in ancestors) {
          corrected_data[nodenum == ancestor, (corrected_column) := TRUE]
        }
      }
    },
    
    "propagate_down" = {
      # If all children are rejected, reject parent
      for (node_id in hierarchy$all_nodes) {
        children <- hierarchy$direct_children[[as.character(node_id)]]
        if (length(children) > 0) {
          child_rejections <- corrected_data[nodenum %in% children, get(corrected_column)]
          if (length(child_rejections) > 0 && all(child_rejections, na.rm = TRUE)) {
            corrected_data[nodenum == node_id, (corrected_column) := TRUE]
          }
        }
      }
    },
    
    "closed_testing" = {
      # Apply proper closed testing (requires implementation of full procedure)
      warning("Full closed testing correction not implemented in this function. Use closed_testing_procedure() instead.")
    },
    
    {
      stop("Unknown correction method: ", correction_method)
    }
  )
  
  return(corrected_data)
}