#' Design Sensitivity Analysis for Hierarchical Testing
#'
#' Implementation of Rosenbaum's design sensitivity analysis to assess
#' robustness of treatment effect findings to unobserved confounding.
#' This extends the standard single-study approach to hierarchical testing.
#'
#' @param idat Individual-level data
#' @param bdat Block-level data
#' @param blockid Block identifier
#' @param formula Test formula
#' @param gamma_range Vector of sensitivity parameters to test
#' @param test_statistic Type of test statistic ("wilcoxon", "sign", "hodges_lehmann")
#' @param alpha Type I error rate
#' @return Design sensitivity analysis results
#' @references
#' Rosenbaum, P. R. (2017). Observation and experiment: an introduction to causal inference.
#' Harvard University Press.
#' @export
design_sensitivity_analysis <- function(idat, bdat, blockid, formula,
                                        gamma_range = seq(1, 3, by = 0.1),
                                        test_statistic = "wilcoxon",
                                        alpha = 0.05) {
  # Parse formula components
  terms <- parse_formula_components(formula)
  outcome <- terms$outcome
  treatment <- terms$treatment

  # Results storage
  sensitivity_results <- data.frame(
    gamma = numeric(0),
    block_id = character(0),
    p_value_lower = numeric(0),
    p_value_upper = numeric(0),
    significant_lower = logical(0),
    significant_upper = logical(0),
    stringsAsFactors = FALSE
  )

  # Analyze each block separately
  unique_blocks <- unique(idat[[blockid]])

  for (block in unique_blocks) {
    block_data <- idat[idat[[blockid]] == block, ]

    # Skip blocks with insufficient data
    if (nrow(block_data) < 4 || length(unique(block_data[[treatment]])) < 2) {
      next
    }

    for (gamma in gamma_range) {
      # Compute sensitivity bounds for this block and gamma
      bounds <- compute_sensitivity_bounds(
        block_data, outcome, treatment, gamma, test_statistic
      )

      sensitivity_results <- rbind(sensitivity_results, data.frame(
        gamma = gamma,
        block_id = block,
        p_value_lower = bounds$p_lower,
        p_value_upper = bounds$p_upper,
        significant_lower = bounds$p_lower <= alpha,
        significant_upper = bounds$p_upper <= alpha,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Compute summary statistics
  summary_stats <- compute_sensitivity_summary(sensitivity_results, alpha)

  return(list(
    sensitivity_results = sensitivity_results,
    summary_stats = summary_stats,
    parameters = list(
      gamma_range = gamma_range,
      test_statistic = test_statistic,
      alpha = alpha
    )
  ))
}

#' Compute sensitivity bounds for a single block
#' @param data Block data
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param gamma Sensitivity parameter
#' @param test_statistic Type of test statistic
#' @return List with upper and lower p-value bounds
compute_sensitivity_bounds <- function(data, outcome, treatment, gamma, test_statistic) {
  y <- data[[outcome]]

  # Handle treatment variable conversion more carefully
  trt_raw <- data[[treatment]]
  if (is.factor(trt_raw)) {
    trt <- as.numeric(trt_raw) - 1
  } else if (is.character(trt_raw)) {
    trt_unique <- unique(trt_raw)
    if (length(trt_unique) == 2) {
      trt <- as.numeric(as.factor(trt_raw)) - 1
    } else {
      trt <- as.numeric(trt_raw)
    }
  } else {
    trt <- as.numeric(trt_raw)
  }

  # Ensure trt is 0/1
  if (length(unique(trt)) == 2 && !all(trt %in% c(0, 1))) {
    trt_min <- min(trt, na.rm = TRUE)
    trt <- trt - trt_min # Make it 0/1
  }

  # Separate treatment and control groups
  y1 <- y[trt == 1]
  y0 <- y[trt == 0]

  n1 <- length(y1)
  n0 <- length(y0)

  if (n1 == 0 || n0 == 0) {
    return(list(p_lower = 1, p_upper = 1))
  }

  # Compute bounds based on test statistic type
  bounds <- switch(test_statistic,
    "wilcoxon" = wilcoxon_sensitivity_bounds(y1, y0, gamma),
    "sign" = sign_test_sensitivity_bounds(y1, y0, gamma),
    "hodges_lehmann" = hodges_lehmann_sensitivity_bounds(y1, y0, gamma),
    stop("Unknown test statistic: ", test_statistic)
  )

  return(bounds)
}

#' Wilcoxon test sensitivity bounds
#' @param y1 Treatment group outcomes
#' @param y0 Control group outcomes
#' @param gamma Sensitivity parameter
#' @return Upper and lower p-value bounds
wilcoxon_sensitivity_bounds <- function(y1, y0, gamma) {
  n1 <- length(y1)
  n0 <- length(y0)

  # Create all pairwise comparisons
  comparisons <- outer(y1, y0, "-")
  w_plus <- sum(comparisons > 0)

  # Total number of comparisons
  total_pairs <- n1 * n0

  # Under no confounding, probability of each comparison being positive is 0.5
  # Under confounding, this ranges from gamma/(1+gamma) to 1/(1+gamma)

  # Compute bounds
  prob_lower <- gamma / (1 + gamma)
  prob_upper <- 1 / (1 + gamma)

  # Use normal approximation for large samples
  if (total_pairs >= 20) {
    mean_lower <- total_pairs * prob_lower
    mean_upper <- total_pairs * prob_upper
    var_est <- total_pairs * 0.25 # Conservative variance estimate

    # Lower bound (most conservative against H0)
    z_lower <- (w_plus - mean_lower) / sqrt(var_est)
    p_lower <- 1 - pnorm(z_lower)

    # Upper bound (least conservative against H0)
    z_upper <- (w_plus - mean_upper) / sqrt(var_est)
    p_upper <- 1 - pnorm(z_upper)
  } else {
    # Use exact computation for small samples (simplified)
    p_lower <- compute_exact_wilcoxon_bound(w_plus, total_pairs, prob_lower)
    p_upper <- compute_exact_wilcoxon_bound(w_plus, total_pairs, prob_upper)
  }

  # Ensure proper ordering: p_lower <= p_upper
  p_lower <- max(0, min(1, p_lower))
  p_upper <- max(0, min(1, p_upper))

  # Fix ordering if necessary
  if (p_lower > p_upper) {
    temp <- p_lower
    p_lower <- p_upper
    p_upper <- temp
  }

  return(list(
    p_lower = p_lower,
    p_upper = p_upper
  ))
}

#' Exact Wilcoxon bound computation (simplified)
#' @param w_plus Number of positive pairwise differences
#' @param total_pairs Total number of pairs
#' @param prob Probability parameter
#' @return P-value bound
compute_exact_wilcoxon_bound <- function(w_plus, total_pairs, prob) {
  # Simplified exact computation using binomial approximation
  p_val <- 1 - pbinom(w_plus - 1, total_pairs, prob)
  return(p_val)
}

#' Sign test sensitivity bounds
#' @param y1 Treatment group outcomes
#' @param y0 Control group outcomes
#' @param gamma Sensitivity parameter
#' @return Upper and lower p-value bounds
sign_test_sensitivity_bounds <- function(y1, y0, gamma) {
  # For sign test, we look at sign of treatment - control within pairs
  # This requires paired data, so we'll approximate with random pairing

  n <- min(length(y1), length(y0))
  if (n == 0) {
    return(list(p_lower = 1, p_upper = 1))
  }

  # Random pairing (in practice, should use actual pairing if available)
  differences <- y1[1:n] - y0[1:n]
  positive_differences <- sum(differences > 0)

  # Sensitivity analysis for sign test
  prob_lower <- gamma / (1 + gamma)
  prob_upper <- 1 / (1 + gamma)

  # Binomial test bounds
  p_lower <- 1 - pbinom(positive_differences - 1, n, prob_lower)
  p_upper <- 1 - pbinom(positive_differences - 1, n, prob_upper)

  return(list(
    p_lower = max(0, min(1, p_lower)),
    p_upper = max(0, min(1, p_upper))
  ))
}

#' Hodges-Lehmann estimator sensitivity bounds
#' @param y1 Treatment group outcomes
#' @param y0 Control group outcomes
#' @param gamma Sensitivity parameter
#' @return Upper and lower p-value bounds
hodges_lehmann_sensitivity_bounds <- function(y1, y0, gamma) {
  # Hodges-Lehmann estimator based on pairwise differences
  all_differences <- outer(y1, y0, "-")
  hl_estimate <- median(as.vector(all_differences))

  # Test whether HL estimator is significantly different from 0
  # This is more complex and requires specialized implementation
  # For now, fall back to Wilcoxon bounds
  return(wilcoxon_sensitivity_bounds(y1, y0, gamma))
}

#' Compute summary statistics for sensitivity analysis
#' @param sensitivity_results Data frame of sensitivity results
#' @param alpha Type I error rate
#' @return Summary statistics
compute_sensitivity_summary <- function(sensitivity_results, alpha) {
  # Find critical gamma values for each block
  block_summaries <- by(sensitivity_results, sensitivity_results$block_id, function(block_data) {
    # Find largest gamma where lower bound is still significant
    sig_lower <- block_data[block_data$significant_lower, ]
    gamma_critical_lower <- if (nrow(sig_lower) > 0) max(sig_lower$gamma) else 1.0

    # Find smallest gamma where upper bound becomes non-significant
    non_sig_upper <- block_data[!block_data$significant_upper, ]
    gamma_critical_upper <- if (nrow(non_sig_upper) > 0) min(non_sig_upper$gamma) else max(block_data$gamma)

    list(
      block_id = unique(block_data$block_id)[1],
      gamma_critical_lower = gamma_critical_lower,
      gamma_critical_upper = gamma_critical_upper,
      robust_to_bias = gamma_critical_lower > 1.0
    )
  })

  # Convert to data frame
  summary_df <- do.call(rbind, lapply(block_summaries, function(x) {
    data.frame(
      block_id = x$block_id,
      gamma_critical_lower = x$gamma_critical_lower,
      gamma_critical_upper = x$gamma_critical_upper,
      robust_to_bias = x$robust_to_bias,
      stringsAsFactors = FALSE
    )
  }))

  # Overall summary
  overall_summary <- list(
    total_blocks = nrow(summary_df),
    blocks_robust_to_bias = sum(summary_df$robust_to_bias),
    mean_gamma_critical = mean(summary_df$gamma_critical_lower),
    median_gamma_critical = median(summary_df$gamma_critical_lower),
    min_gamma_critical = min(summary_df$gamma_critical_lower),
    max_gamma_critical = max(summary_df$gamma_critical_lower)
  )

  return(list(
    block_summaries = summary_df,
    overall_summary = overall_summary
  ))
}

#' Plot sensitivity analysis results
#' @param sensitivity_analysis Results from design_sensitivity_analysis
#' @param blocks Vector of block IDs to plot (NULL for all)
#' @return ggplot object
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_hline facet_wrap labs theme_minimal theme element_text
#' @export
plot_design_sensitivity <- function(sensitivity_analysis, blocks = NULL) {
  results <- sensitivity_analysis$sensitivity_results

  if (!is.null(blocks)) {
    results <- results[results$block_id %in% blocks, ]
  }

  # Create the plot
  p <- ggplot(results, aes(x = gamma)) +
    geom_ribbon(aes(ymin = p_value_lower, ymax = p_value_upper, fill = block_id),
      alpha = 0.3
    ) +
    geom_hline(
      yintercept = sensitivity_analysis$parameters$alpha,
      linetype = "dashed", color = "red"
    ) +
    facet_wrap(~block_id, scales = "free_y") +
    labs(
      x = "Sensitivity Parameter (Gamma)",
      y = "P-value Bounds",
      title = "Design Sensitivity Analysis",
      subtitle = paste("Red line indicates alpha =", sensitivity_analysis$parameters$alpha)
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  return(p)
}

#' Parse formula components
#' @param formula Formula object
#' @return List with parsed components
parse_formula_components <- function(formula) {
  formula_str <- deparse(formula)
  parts <- strsplit(formula_str, "~")[[1]]

  outcome <- trimws(parts[1])
  rhs <- trimws(parts[2])

  # Handle blocking structure
  if (grepl("\\|", rhs)) {
    treatment_part <- trimws(strsplit(rhs, "\\|")[[1]][1])
    block_part <- trimws(strsplit(rhs, "\\|")[[1]][2])
  } else {
    treatment_part <- rhs
    block_part <- NULL
  }

  return(list(
    outcome = outcome,
    treatment = treatment_part,
    block = block_part
  ))
}
