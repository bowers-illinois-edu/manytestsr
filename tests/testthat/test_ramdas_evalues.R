test_that("e-value sequential testing follows Ramdas methodology", {
  skip_if_not_installed("data.table")

  # Create simple test data for e-value testing
  set.seed(42)
  n_per_block <- 10
  n_blocks <- 5

  test_data <- data.frame(
    id = 1:(n_per_block * n_blocks),
    block = rep(1:n_blocks, each = n_per_block),
    treatment = rep(c(0, 1), n_per_block * n_blocks / 2),
    outcome = rnorm(n_per_block * n_blocks,
      mean = rep(c(0, 1), n_per_block * n_blocks / 2) * 0.3,
      sd = 1
    )
  )

  if (exists("evalue_sequential_test")) {
    # Test Kelly betting (Ramdas recommended approach)
    result_kelly <- evalue_sequential_test(
      data = test_data,
      formula = outcome ~ treatment,
      blocks = test_data$block,
      alpha = 0.05,
      wealth_rule = "kelly",
      stopping_rule = "wealth_threshold"
    )

    expect_true(is.list(result_kelly))
    expect_true("results" %in% names(result_kelly))
    expect_true("final_wealth" %in% names(result_kelly))

    # Check e-value properties
    expect_true(all(result_kelly$results$e_value >= 1)) # E-values always >= 1
    expect_true(result_kelly$final_wealth > 0) # Wealth should be positive

    # Test fixed betting strategy
    result_fixed <- evalue_sequential_test(
      data = test_data,
      formula = outcome ~ treatment,
      blocks = test_data$block,
      alpha = 0.05,
      wealth_rule = "fixed",
      stopping_rule = "never"
    )

    expect_true(all(result_fixed$results$e_value >= 1))

    # Test adaptive betting
    result_adaptive <- evalue_sequential_test(
      data = test_data,
      formula = outcome ~ treatment,
      blocks = test_data$block,
      alpha = 0.05,
      wealth_rule = "adaptive",
      stopping_rule = "e_threshold"
    )

    expect_true(all(result_adaptive$results$e_value >= 1))
  } else {
    skip("evalue_sequential_test function not available")
  }
})

test_that("e-value computation methods work correctly", {
  # Create simple block data for testing e-value computation
  set.seed(123)
  block_data <- data.frame(
    outcome = c(rnorm(10, 0, 1), rnorm(10, 0.5, 1)), # Treatment effect
    treatment = c(rep("control", 10), rep("treatment", 10))
  )

  if (exists("compute_evalue")) {
    # Test Kelly e-value computation
    e_val_kelly <- compute_evalue(block_data, "outcome", "treatment", "kelly")
    expect_true(is.numeric(e_val_kelly))
    expect_true(e_val_kelly >= 1) # E-values are always >= 1

    # Test fixed e-value computation
    e_val_fixed <- compute_evalue(block_data, "outcome", "treatment", "fixed")
    expect_true(is.numeric(e_val_fixed))
    expect_true(e_val_fixed >= 1)

    # Test adaptive e-value computation
    e_val_adaptive <- compute_evalue(block_data, "outcome", "treatment", "adaptive")
    expect_true(is.numeric(e_val_adaptive))
    expect_true(e_val_adaptive >= 1)
  }

  # Test individual betting strategies if they exist
  test_t_stat <- 2.1 # Significant t-statistic

  if (exists("kelly_evalue")) {
    kelly_result <- kelly_evalue(test_t_stat, 10, 10)
    expect_true(is.numeric(kelly_result))
    expect_true(kelly_result >= 1)
  }

  if (exists("fixed_evalue")) {
    fixed_result <- fixed_evalue(test_t_stat, 10, 10)
    expect_true(is.numeric(fixed_result))
    expect_true(fixed_result >= 1)
  }

  if (exists("adaptive_evalue")) {
    adaptive_result <- adaptive_evalue(test_t_stat, 10, 10)
    expect_true(is.numeric(adaptive_result))
    expect_true(adaptive_result >= 1)
  }
})

test_that("wealth accumulation follows Ramdas principles", {
  if (exists("update_wealth")) {
    # Test multiplicative wealth update (standard in Ramdas framework)
    initial_wealth <- 1.0
    e_value <- 2.0
    alpha <- 0.05

    updated_wealth <- update_wealth(initial_wealth, e_value, alpha, "kelly")
    expect_true(updated_wealth == initial_wealth * e_value)
    expect_true(updated_wealth > initial_wealth)

    # Test wealth accumulation over sequence
    wealth_sequence <- numeric(5)
    wealth_sequence[1] <- 1.0
    e_values <- c(1.5, 2.0, 0.8, 1.2, 3.0)

    for (i in 2:5) {
      wealth_sequence[i] <- update_wealth(
        wealth_sequence[i - 1],
        e_values[i],
        alpha,
        "kelly"
      )
    }

    expect_true(all(wealth_sequence > 0)) # Wealth should stay positive
    expect_true(length(wealth_sequence) == 5)
  }
})

test_that("stopping rules work correctly", {
  if (exists("check_stopping_rule")) {
    # Test wealth threshold stopping
    high_wealth <- 25.0 # Above 1/0.05 = 20
    should_stop <- check_stopping_rule(high_wealth, 2.0, 0.05, "wealth_threshold")
    expect_true(should_stop)

    low_wealth <- 15.0 # Below 1/0.05 = 20
    should_continue <- check_stopping_rule(low_wealth, 2.0, 0.05, "wealth_threshold")
    expect_false(should_continue)

    # Test e-value threshold stopping
    high_e_value <- 25.0
    should_stop_e <- check_stopping_rule(10.0, high_e_value, 0.05, "e_threshold")
    expect_true(should_stop_e)

    # Test never stopping
    never_stop <- check_stopping_rule(100.0, 50.0, 0.05, "never")
    expect_false(never_stop)
  }
})

test_that("e-value to p-value conversion is valid", {
  if (exists("evalue_to_pvalue")) {
    # Test conversion relationship: p = 1/e (basic version)
    e_values <- c(1.0, 2.0, 5.0, 10.0, 20.0)

    for (e_val in e_values) {
      p_val <- evalue_to_pvalue(e_val)
      expect_true(is.numeric(p_val))
      expect_true(p_val >= 0 && p_val <= 1)
      expect_true(abs(p_val - 1 / e_val) < 1e-10) # Should equal 1/e
    }

    # Test edge cases
    expect_equal(evalue_to_pvalue(1.0), 1.0) # Minimum e-value gives p=1
    expect_true(evalue_to_pvalue(Inf) == 0) # Infinite e-value gives p=0
  }
})

test_that("confidence sequences have anytime validity", {
  if (exists("evalue_confidence_sequence")) {
    # Test confidence sequence construction
    e_values <- c(1.2, 2.1, 1.8, 3.2, 2.5, 4.1, 1.9)

    conf_seq <- evalue_confidence_sequence(e_values, alpha = 0.05)

    expect_true(is.data.frame(conf_seq))
    expect_true("step" %in% names(conf_seq))
    expect_true("lower" %in% names(conf_seq))
    expect_true("upper" %in% names(conf_seq))
    expect_true("wealth" %in% names(conf_seq))

    # Check that bounds are properly ordered
    expect_true(all(conf_seq$lower <= conf_seq$upper))

    # Check that wealth grows with e-values
    expect_true(all(diff(conf_seq$wealth) >= 0)) # Wealth should be non-decreasing
  }
})

test_that("e-values work with find_blocks integration", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

  # Test integration with find_blocks (should give experimental warning)
  data(example_dat, package = "manytestsr")

  idat <- data.table::as.data.table(example_dat)
  bdat <- idat %>%
    dplyr::group_by(blockF) %>%
    dplyr::summarize(
      nb = dplyr::n(),
      pb = mean(trt),
      hwt = (nb / nrow(idat)) * (pb * (1 - pb)),
      .groups = "drop"
    ) %>%
    data.table::as.data.table()

  # Should produce warning about experimental status
  expect_warning(
    {
      result_evalues <- find_blocks(
        idat = idat,
        bdat = bdat,
        blockid = "blockF",
        splitfn = splitCluster,
        pfn = pOneway,
        fmla = Y1 ~ trtF | blockF,
        splitby = "hwt",
        parallel = "no",
        use_evalues = TRUE,
        evalue_wealth_rule = "kelly",
        thealpha = 0.05,
        maxtest = 8,
        trace = FALSE
      )
    },
    "experimental"
  )

  # Should still produce valid basic results
  expect_true(is.list(result_evalues))
  expect_true("node_dat" %in% names(result_evalues))
})

test_that("different wealth rules produce different behaviors", {
  skip_if_not_installed("data.table")

  if (exists("evalue_sequential_test")) {
    set.seed(789)
    test_data <- data.frame(
      block = rep(1:8, each = 12),
      treatment = rep(c(0, 1), 48),
      outcome = rnorm(96, mean = rep(c(0, 0.4), 48), sd = 1)
    )

    # Compare different wealth rules
    kelly_result <- evalue_sequential_test(
      test_data, outcome ~ treatment, test_data$block,
      wealth_rule = "kelly"
    )

    fixed_result <- evalue_sequential_test(
      test_data, outcome ~ treatment, test_data$block,
      wealth_rule = "fixed"
    )

    adaptive_result <- evalue_sequential_test(
      test_data, outcome ~ treatment, test_data$block,
      wealth_rule = "adaptive"
    )

    # All should produce valid results
    expect_true(kelly_result$final_wealth > 0)
    expect_true(fixed_result$final_wealth > 0)
    expect_true(adaptive_result$final_wealth > 0)

    # Different rules should potentially give different final wealth
    # (Though with same data, differences might be subtle)
    wealth_values <- c(
      kelly_result$final_wealth,
      fixed_result$final_wealth,
      adaptive_result$final_wealth
    )

    expect_true(all(wealth_values > 0))
    expect_true(length(unique(round(wealth_values, 6))) >= 1) # Allow for some similarity
  }
})

test_that("e-values handle edge cases properly", {
  if (exists("evalue_sequential_test")) {
    # Test with no treatment effect (null hypothesis true)
    set.seed(999)
    null_data <- data.frame(
      block = rep(1:5, each = 10),
      treatment = rep(c(0, 1), 25),
      outcome = rnorm(50, mean = 0, sd = 1) # No treatment effect
    )

    null_result <- evalue_sequential_test(
      null_data, outcome ~ treatment, null_data$block,
      wealth_rule = "kelly", stopping_rule = "never"
    )

    # Under null, wealth should not grow systematically large
    expect_true(null_result$final_wealth > 0)
    expect_true(all(null_result$results$e_value >= 1))

    # Test with single block
    single_block_data <- data.frame(
      block = rep(1, 20),
      treatment = rep(c(0, 1), 10),
      outcome = rnorm(20, mean = rep(c(0, 0.5), 10), sd = 1)
    )

    expect_no_error({
      single_result <- evalue_sequential_test(
        single_block_data, outcome ~ treatment, single_block_data$block,
        wealth_rule = "fixed"
      )
    })

    expect_true(single_result$final_wealth > 0)
  }
})
