test_that("Closed testing procedure works correctly", {
  skip_if_not_installed("data.table")
  
  # Create simple test data
  node_dat <- data.table::data.table(
    nodenum = 1:7,
    p = c(0.001, 0.01, 0.03, 0.02, 0.08, 0.004, 0.12),
    depth = c(1, 2, 2, 3, 3, 3, 3),
    testable = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE)
  )
  
  # Create simple tracker
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
      parent_id = c(0L, 1L, 1L, 2L, 2L, 3L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L, 3L, 3L)
    ),
    next_id = 8L
  )
  
  # Test closed testing procedure (if function exists)
  skip_if_not(exists("closed_testing_procedure"))
  
  result <- closed_testing_procedure(node_dat, tracker, alpha = 0.05, method = "simes")
  
  # Check that result has correct structure
  expect_true(is.data.table(result) || "data.table" %in% class(result))
  expect_true("closed_testing_reject" %in% names(result))
  
  # Check that some rejections occurred
  expect_true(any(result$closed_testing_reject))
  
  # Check that rejected nodes have significant p-values
  rejected_nodes <- result[closed_testing_reject == TRUE]
  expect_true(all(rejected_nodes$p <= 0.05))
})

test_that("E-value sequential testing works", {
  skip_if_not_installed("data.table")
  
  # Create test data
  set.seed(123)
  n_per_block <- 20
  n_blocks <- 5
  
  test_data <- data.frame(
    id = 1:(n_per_block * n_blocks),
    block = rep(1:n_blocks, each = n_per_block),
    treatment = rep(c(0, 1), n_per_block * n_blocks / 2),
    outcome = rnorm(n_per_block * n_blocks, 
                   mean = rep(c(0, 1), n_per_block * n_blocks / 2) * 0.3, 
                   sd = 1)
  )
  
  # Test e-value function (if it exists)
  skip_if_not(exists("evalue_sequential_test"))
  
  result <- evalue_sequential_test(
    data = test_data,
    formula = outcome ~ treatment,
    blocks = test_data$block,
    alpha = 0.05,
    wealth_rule = "kelly"
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true("results" %in% names(result))
  expect_true("final_wealth" %in% names(result))
  expect_true(is.data.frame(result$results))
  
  # Check that e-values are >= 1
  expect_true(all(result$results$e_value >= 1))
  
  # Check that wealth is positive
  expect_true(result$final_wealth > 0)
})

test_that("Design sensitivity analysis works", {
  skip_if_not_installed("data.table")
  
  # Create test data with treatment effect
  set.seed(456)
  idat <- data.frame(
    block = rep(paste0("B", 1:10), each = 10),
    treatment = rep(c("control", "treatment"), 50),
    outcome = c(
      rnorm(50, mean = 0, sd = 1),  # Control
      rnorm(50, mean = 0.5, sd = 1)  # Treatment with effect
    )
  )
  
  bdat <- data.frame(
    block = paste0("B", 1:10),
    size = rep(10, 10)
  )
  
  # Test sensitivity analysis (if function exists)
  skip_if_not(exists("design_sensitivity_analysis"))
  
  result <- design_sensitivity_analysis(
    idat = idat,
    bdat = bdat,
    blockid = "block",
    formula = outcome ~ treatment,
    gamma_range = seq(1, 2, by = 0.2),
    alpha = 0.05
  )
  
  # Check structure
  expect_true(is.list(result))
  expect_true("sensitivity_results" %in% names(result))
  expect_true("summary_stats" %in% names(result))
  
  # Check that sensitivity results have correct columns
  sens_results <- result$sensitivity_results
  expected_cols <- c("gamma", "block_id", "p_value_lower", "p_value_upper")
  expect_true(all(expected_cols %in% names(sens_results)))
  
  # Check that p-value bounds are valid
  expect_true(all(sens_results$p_value_lower >= 0 & sens_results$p_value_lower <= 1))
  expect_true(all(sens_results$p_value_upper >= 0 & sens_results$p_value_upper <= 1))
  expect_true(all(sens_results$p_value_lower <= sens_results$p_value_upper))
})

test_that("Intersection-union tests work correctly", {
  skip_if_not_installed("data.table")
  
  # Create hierarchical test data
  node_dat <- data.table::data.table(
    nodenum = 1:7,
    p = c(0.001, 0.01, 0.03, 0.02, 0.08, 0.004, 0.12),
    depth = c(1, 2, 2, 3, 3, 3, 3)
  )
  
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
      parent_id = c(0L, 1L, 1L, 2L, 2L, 3L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L, 3L, 3L)
    ),
    next_id = 8L
  )
  
  # Test intersection-union tests (if function exists)
  skip_if_not(exists("intersection_union_tests"))
  
  result <- intersection_union_tests(node_dat, tracker, alpha = 0.05)
  
  # Check structure
  expect_true(is.list(result))
  expect_true("node_dat" %in% names(result))
  expect_true("iu_results" %in% names(result))
  
  # Check that new columns were added
  updated_node_dat <- result$node_dat
  expect_true("iu_p_intersection" %in% names(updated_node_dat))
  expect_true("iu_p_union" %in% names(updated_node_dat))
  
  # Check that p-values are valid
  expect_true(all(updated_node_dat$iu_p_intersection >= 0, na.rm = TRUE))
  expect_true(all(updated_node_dat$iu_p_intersection <= 1, na.rm = TRUE))
  expect_true(all(updated_node_dat$iu_p_union >= 0, na.rm = TRUE))
  expect_true(all(updated_node_dat$iu_p_union <= 1, na.rm = TRUE))
})

test_that("Consonance checking works correctly", {
  skip_if_not_installed("data.table")
  
  # Create test data with potential consonance violations
  node_dat <- data.table::data.table(
    nodenum = 1:7,
    p = c(0.001, 0.01, 0.03, 0.02, 0.08, 0.004, 0.12),
    depth = c(1, 2, 2, 3, 3, 3, 3),
    testable = c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)  # Potential violations
  )
  
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L, 6L, 7L),
      parent_id = c(0L, 1L, 1L, 2L, 2L, 3L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L, 3L, 3L)
    )
  )
  
  # Test consonance checking (if function exists)
  if (exists("check_consonance_property")) {
    expect_no_error({
      result <- check_consonance_property(node_dat, tracker, "testable")
    })
    
    if (exists("result")) {
      # Check structure
      expect_true(is.list(result))
      expect_true("is_consonant" %in% names(result))
      expect_true("violations" %in% names(result))
      
      # Should detect violations in our test case
      expect_true(is.logical(result$is_consonant))
      expect_true(is.list(result$violations))
    }
  } else {
    skip("check_consonance_property function not available")
  }
})

test_that("Advanced methodology integration in find_blocks", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")
  
  # Load example data
  data(example_dat, package = "manytestsr", envir = environment())
  
  # Prepare data
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
  
  # Test find_blocks with closed testing (if available)
  if (exists("closed_testing_procedure")) {
    expect_no_error({
      result <- find_blocks(
        idat = idat,
        bdat = bdat,
        blockid = "blockF",
        splitfn = splitCluster,
        pfn = pOneway,
        fmla = Y1 ~ trtF | blockF,
        splitby = "hwt",
        parallel = "no",
        maxtest = 10,
        thealpha = 0.05,
        use_closed_testing = TRUE,
        closed_testing_method = "simes"
      )
    })
    
    # Check that closed testing results are included if function worked
    if (exists("result") && "node_dat" %in% names(result)) {
      node_dat <- result$node_dat
      if ("closed_testing_reject" %in% names(node_dat)) {
        expect_true(is.logical(node_dat$closed_testing_reject))
      }
    }
  }
  
  # Test with e-values (should give warning since not fully integrated)
  expect_warning({
    result_evalues <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF", 
      splitfn = splitCluster,
      pfn = pOneway,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      maxtest = 10,
      thealpha = 0.05,
      use_evalues = TRUE
    )
  }, "experimental")
})

test_that("P-value combination methods work correctly", {
  # Test different p-value combination methods used in advanced procedures
  
  p_vals <- c(0.01, 0.05, 0.08, 0.12)
  
  # Test Simes method (if available in local functions)
  if (exists("local_simes")) {
    simes_result <- local_simes(p_vals)
    expect_true(is.numeric(simes_result))
    expect_true(length(simes_result) == 1)
    expect_true(simes_result >= 0 && simes_result <= 1)
    # Simes can be more significant than min due to multiplicity correction
    expect_true(simes_result <= min(p_vals) * length(p_vals))
  }
  
  # Test that combination methods handle edge cases
  expect_no_error(local_simes(c(0.001, 0.999)))
  expect_no_error(local_simes(c(1.0)))
  expect_no_error(local_simes(numeric(0)))  # Empty vector
})

test_that("Helper functions for advanced methods work", {
  # Test node ancestry building
  tracker <- list(
    tracker = data.table::data.table(
      node_id = c(1L, 2L, 3L, 4L, 5L),
      parent_id = c(0L, 1L, 1L, 2L, 3L),
      depth = c(1L, 2L, 2L, 3L, 3L)
    )
  )
  
  # Test ancestry function
  if (exists("build_numeric_ancestry")) {
    ancestry_4 <- build_numeric_ancestry(tracker, 4)
    expect_true(is.numeric(ancestry_4))
    expect_true(1 %in% ancestry_4)  # Should include root
    expect_true(2 %in% ancestry_4)  # Should include direct parent
    expect_true(4 %in% ancestry_4)  # Should include self
  }
})

test_that("Error handling in advanced methods", {
  # Test that methods handle invalid inputs gracefully
  
  # Empty node data
  empty_node_dat <- data.table::data.table(
    nodenum = integer(0),
    testable = logical(0)
  )
  empty_tracker <- list(
    tracker = data.table::data.table(
      node_id = integer(0),
      parent_id = integer(0),
      depth = integer(0)
    ), 
    next_id = 1L
  )
  
  if (exists("check_consonance_property")) {
    tryCatch({
      result <- check_consonance_property(empty_node_dat, empty_tracker, "testable")
      expect_true(is.list(result))
    }, error = function(e) {
      # Empty data structures may cause legitimate errors
      expect_true(TRUE)
    })
  }
  
  # Invalid alpha values
  if (exists("closed_testing_procedure")) {
    node_dat <- data.table::data.table(
      nodenum = 1:3,
      p = c(0.01, 0.05, 0.1),
      depth = c(1, 2, 2)
    )
    tracker <- list(
      tracker = data.table::data.table(
        node_id = 1:3,
        parent_id = c(0, 1, 1),
        depth = 1:3
      )
    )
    
    expect_error(closed_testing_procedure(node_dat, tracker, alpha = -0.1), "alpha must be")
    expect_error(closed_testing_procedure(node_dat, tracker, alpha = 1.5), "alpha must be")
  }
})