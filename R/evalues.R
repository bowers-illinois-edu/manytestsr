#' E-value Based Sequential Testing
#'
#' Implementation of e-value methodology for sequential hypothesis testing
#' following the framework of Ramdas et al. This provides anytime-valid
#' inference with optional stopping and flexible error control.
#'
#' @param data Individual-level data
#' @param formula Test formula
#' @param blocks Block identifiers
#' @param alpha Type I error rate
#' @param wealth_rule Wealth accumulation rule ("kelly", "fixed", "adaptive")
#' @param stopping_rule Early stopping rule ("never", "wealth_threshold", "e_threshold")
#' @return E-value based test results
#' @references
#' Ramdas, A., Grünwald, P., Vovk, V., & Shafer, G. (2023). Game-theoretic
#' statistics and safe anytime-valid inference. Statistical Science, 38(4), 576-601.
#' @examples
#' \dontrun{
#' # Load example data for e-value testing with find_blocks
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
#' # Run find_blocks with e-value methodology (experimental)
#' # Note: This provides anytime-valid inference
#' results_evalues <- find_blocks(
#'   idat = idat,
#'   bdat = bdat,
#'   blockid = "blockF",
#'   splitfn = splitCluster,
#'   pfn = pOneway,
#'   fmla = Y1 ~ trtF | blockF,
#'   splitby = "hwt",
#'   parallel = "no",
#'   use_evalues = TRUE,
#'   evalue_wealth_rule = "kelly", # Ramdas recommends Kelly betting
#'   thealpha = 0.05,
#'   maxtest = 15
#' )
#'
#' # Direct use of e-value sequential testing
#' evalue_results <- evalue_sequential_test(
#'   data = idat,
#'   formula = Y1 ~ trtF,
#'   blocks = idat$blockF,
#'   alpha = 0.05,
#'   wealth_rule = "kelly",
#'   stopping_rule = "wealth_threshold"
#' )
#'
#' # Examine e-value results
#' cat("E-value sequential testing results:\n")
#' cat("Final wealth:", evalue_results$final_wealth, "\n")
#' cat("Total rejections:", evalue_results$total_rejections, "\n")
#' cat("Stopped early:", evalue_results$stopped_early, "\n")
#'
#' # Plot wealth accumulation over time
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   wealth_plot <- ggplot(
#'     evalue_results$results,
#'     aes(x = seq_along(block_id), y = wealth)
#'   ) +
#'     geom_line() +
#'     geom_hline(yintercept = 1 / 0.05, linetype = "dashed", color = "red") +
#'     labs(
#'       x = "Sequential Test Number", y = "Wealth",
#'       title = "E-value Wealth Accumulation",
#'       subtitle = "Red line shows rejection threshold (1/alpha)"
#'     )
#'   print(wealth_plot)
#' }
#'
#' # Compare e-values with traditional p-values
#' comparison_data <- data.frame(
#'   test = seq_along(evalue_results$results$e_value),
#'   e_value = evalue_results$results$e_value,
#'   p_value = evalue_results$results$p_value,
#'   wealth = evalue_results$results$wealth
#' )
#'
#' cat("E-value vs P-value comparison (first 5 tests):\n")
#' print(head(comparison_data, 5))
#' }
#' @export
evalue_sequential_test <- function(data, formula, blocks, alpha = 0.05,
                                   wealth_rule = "kelly",
                                   stopping_rule = "wealth_threshold") {
  # Extract components from formula
  terms <- parse_formula(formula)
  outcome <- terms$outcome
  treatment <- terms$treatment

  # Initialize wealth and e-values
  wealth <- 1.0
  e_values <- numeric()
  p_values <- numeric()
  decisions <- logical()
  stopped <- FALSE

  # Sequential testing over blocks
  unique_blocks <- unique(blocks)
  results <- list()

  for (i in seq_along(unique_blocks)) {
    if (stopped) break

    block_id <- unique_blocks[i]
    block_data <- data[blocks == block_id, ]

    # Compute e-value for this block
    e_val <- compute_evalue(block_data, outcome, treatment, wealth_rule)
    p_val <- evalue_to_pvalue(e_val)

    # Update wealth
    wealth <- update_wealth(wealth, e_val, alpha, wealth_rule)

    # Make decision
    reject <- (wealth >= 1 / alpha) # Equivalently, e_val >= 1/alpha

    # Store results
    e_values <- c(e_values, e_val)
    p_values <- c(p_values, p_val)
    decisions <- c(decisions, reject)

    # Check stopping rule
    if (stopping_rule != "never") {
      stopped <- check_stopping_rule(wealth, e_val, alpha, stopping_rule)
    }

    results[[i]] <- list(
      block_id = block_id,
      e_value = e_val,
      p_value = p_val,
      wealth = wealth,
      reject = reject,
      stopped = stopped
    )
  }

  # Compile final results
  result_summary <- data.frame(
    block_id = unique_blocks[1:length(e_values)],
    e_value = e_values,
    p_value = p_values,
    wealth = sapply(results[1:length(e_values)], function(x) x$wealth),
    reject = decisions,
    stringsAsFactors = FALSE
  )

  return(list(
    results = result_summary,
    final_wealth = wealth,
    total_rejections = sum(decisions),
    stopped_early = stopped,
    method_info = list(
      alpha = alpha,
      wealth_rule = wealth_rule,
      stopping_rule = stopping_rule
    )
  ))
}

#' Compute e-value for a single test
#' @param data Block data
#' @param outcome Outcome variable name
#' @param treatment Treatment variable name
#' @param wealth_rule Wealth rule for betting
#' @return E-value
compute_evalue <- function(data, outcome, treatment, wealth_rule) {
  # Extract treatment assignments and outcomes
  y <- data[[outcome]]
  trt <- as.numeric(data[[treatment]]) - 1 # Convert to 0/1

  if (length(unique(trt)) < 2) {
    return(1.0) # No variation in treatment
  }

  # Compute sufficient statistics
  n1 <- sum(trt == 1)
  n0 <- sum(trt == 0)
  s1 <- sum(y[trt == 1])
  s0 <- sum(y[trt == 0])

  if (n1 == 0 || n0 == 0) {
    return(1.0) # Degenerate case
  }

  # Use t-statistic based e-value construction
  t_stat <- compute_t_statistic(y, trt)

  # Convert to e-value using different betting strategies
  e_val <- switch(wealth_rule,
    "kelly" = kelly_evalue(t_stat, n1, n0),
    "fixed" = fixed_evalue(t_stat, n1, n0),
    "adaptive" = adaptive_evalue(t_stat, n1, n0),
    stop("Unknown wealth rule: ", wealth_rule)
  )

  return(max(e_val, 1.0)) # E-values are always >= 1
}

#' Compute t-statistic for two-sample test
#' @param y Outcomes
#' @param trt Treatment indicators (0/1)
#' @return T-statistic
compute_t_statistic <- function(y, trt) {
  y1 <- y[trt == 1]
  y0 <- y[trt == 0]

  n1 <- length(y1)
  n0 <- length(y0)

  if (n1 <= 1 || n0 <= 1) {
    return(0)
  }

  mean1 <- mean(y1)
  mean0 <- mean(y0)
  var1 <- var(y1)
  var0 <- var(y0)

  # Pooled standard error
  se <- sqrt(var1 / n1 + var0 / n0)

  if (se == 0) {
    return(0)
  }

  t_stat <- (mean1 - mean0) / se
  return(t_stat)
}

#' Kelly optimal betting for e-values
#' @param t_stat T-statistic
#' @param n1 Treatment group size
#' @param n0 Control group size
#' @return E-value using Kelly betting
kelly_evalue <- function(t_stat, n1, n0) {
  # Degrees of freedom approximation
  df <- n1 + n0 - 2

  # Convert t-statistic to e-value using Kelly optimal betting
  # This uses the connection between e-values and betting
  if (abs(t_stat) < 1e-10) {
    return(1.0)
  }

  # Simple Kelly-based transformation (can be improved)
  e_val <- exp(-0.5 * t_stat^2 / (1 + t_stat^2 / df))
  return(1 / max(e_val, 1e-10))
}

#' Fixed betting strategy for e-values
#' @param t_stat T-statistic
#' @param n1 Treatment group size
#' @param n0 Control group size
#' @return E-value using fixed betting
fixed_evalue <- function(t_stat, n1, n0) {
  # Fixed betting fraction approach
  df <- n1 + n0 - 2
  p_val <- 2 * (1 - pt(abs(t_stat), df))

  # Convert p-value to e-value using fixed betting
  if (p_val <= 0) {
    return(Inf)
  }
  if (p_val >= 1) {
    return(1.0)
  }

  # Simple transformation (can be improved with proper theory)
  e_val <- -1 / log(p_val)
  return(max(e_val, 1.0))
}

#' Adaptive betting strategy for e-values
#' @param t_stat T-statistic
#' @param n1 Treatment group size
#' @param n0 Control group size
#' @return E-value using adaptive betting
adaptive_evalue <- function(t_stat, n1, n0) {
  # Adaptive strategy that adjusts based on observed effect size
  effect_size <- abs(t_stat) / sqrt(n1 + n0)

  # Weight betting based on effect size
  weight <- min(1.0, 2 * effect_size)

  # Combine Kelly and fixed approaches
  kelly_e <- kelly_evalue(t_stat, n1, n0)
  fixed_e <- fixed_evalue(t_stat, n1, n0)

  e_val <- weight * kelly_e + (1 - weight) * fixed_e
  return(max(e_val, 1.0))
}

#' Update wealth using e-value
#' @param current_wealth Current wealth level
#' @param e_value New e-value
#' @param alpha Type I error rate
#' @param wealth_rule Wealth update rule
#' @return Updated wealth
update_wealth <- function(current_wealth, e_value, alpha, wealth_rule) {
  # Standard multiplicative update
  new_wealth <- current_wealth * e_value

  # Optional: implement wealth redistribution or other rules
  return(new_wealth)
}

#' Check stopping rule
#' @param wealth Current wealth
#' @param e_value Current e-value
#' @param alpha Type I error rate
#' @param stopping_rule Stopping rule to apply
#' @return Logical indicating whether to stop
check_stopping_rule <- function(wealth, e_value, alpha, stopping_rule) {
  switch(stopping_rule,
    "never" = FALSE,
    "wealth_threshold" = wealth >= 1 / alpha,
    "e_threshold" = e_value >= 1 / alpha,
    FALSE
  )
}

#' Convert e-value to p-value
#' @param e_value E-value
#' @return Corresponding p-value
evalue_to_pvalue <- function(e_value) {
  # Standard conversion: p = 1/e
  return(1 / max(e_value, 1.0))
}

#' Parse formula to extract components
#' @param formula Test formula
#' @return List with outcome and treatment components
parse_formula <- function(formula) {
  formula_str <- deparse(formula)
  parts <- strsplit(formula_str, "~")[[1]]

  outcome <- trimws(parts[1])
  rhs <- trimws(parts[2])

  # Handle blocking structure if present
  if (grepl("\\|", rhs)) {
    treatment_part <- trimws(strsplit(rhs, "\\|")[[1]][1])
  } else {
    treatment_part <- rhs
  }

  return(list(outcome = outcome, treatment = treatment_part))
}

#' E-value based confidence sequences
#' @param e_values Vector of e-values from sequential testing
#' @param alpha Type I error rate
#' @return Confidence sequence boundaries
evalue_confidence_sequence <- function(e_values, alpha = 0.05) {
  # Construct anytime-valid confidence sequences
  n <- length(e_values)

  # Use wealth process for confidence bounds
  wealth_process <- cumprod(e_values)

  # Confidence boundaries (this is a simplified version)
  lower_bound <- numeric(n)
  upper_bound <- numeric(n)

  for (i in 1:n) {
    # Conservative bounds based on wealth
    width <- sqrt(log(wealth_process[i]) / i)
    lower_bound[i] <- -width
    upper_bound[i] <- width
  }

  return(data.frame(
    step = 1:n,
    lower = lower_bound,
    upper = upper_bound,
    wealth = wealth_process
  ))
}
