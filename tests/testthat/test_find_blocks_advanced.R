test_that("find_blocks works with closed testing procedure", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

  # Load and prepare example data
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

  # Test find_blocks with closed testing
  expect_no_error({
    results_closed <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pOneway,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_closed_testing = TRUE,
      closed_testing_method = "simes",
      thealpha = 0.05,
      maxtest = 10,
      trace = FALSE
    )
  })

  # Check basic structure
  expect_true(is.list(results_closed))
  expect_true("bdat" %in% names(results_closed))
  expect_true("node_dat" %in% names(results_closed))

  # Check if closed testing was applied (depends on function existence)
  if ("closed_testing_reject" %in% names(results_closed$node_dat)) {
    expect_true(is.logical(results_closed$node_dat$closed_testing_reject))

    # Check that rejected nodes have significant p-values
    rejected_nodes <- results_closed$node_dat[closed_testing_reject == TRUE]
    if (nrow(rejected_nodes) > 0) {
      expect_true(all(rejected_nodes$p <= 0.05, na.rm = TRUE))
    }
  }
})

test_that("find_blocks works with different closed testing methods", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  methods_to_test <- c("simes", "fisher", "min")

  for (method in methods_to_test) {
    expect_no_error({
      result <- find_blocks(
        idat = idat,
        bdat = bdat,
        blockid = "blockF",
        splitfn = splitCluster,
        pfn = pIndepDist,
        fmla = Y1 ~ trtF | blockF,
        splitby = "hwt",
        parallel = "no",
        use_closed_testing = TRUE,
        closed_testing_method = method,
        thealpha = 0.05,
        maxtest = 8,
        trace = FALSE
      )
    })

    # Basic checks
    expect_true(is.list(result))
    expect_true("node_dat" %in% names(result))
  }
})

test_that("find_blocks gives warning for e-values (experimental feature)", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  # Should give warning since e-value integration is experimental
  expect_warning(
    {
      results_evalues <- find_blocks(
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

  # Should still produce valid results
  expect_true(is.list(results_evalues))
  expect_true("bdat" %in% names(results_evalues))
  expect_true("node_dat" %in% names(results_evalues))
})

test_that("closed testing vs traditional comparison works", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  # Traditional approach
  results_traditional <- find_blocks(
    idat = idat,
    bdat = bdat,
    blockid = "blockF",
    splitfn = splitCluster,
    pfn = pOneway,
    fmla = Y1 ~ trtF | blockF,
    splitby = "hwt",
    parallel = "no",
    use_closed_testing = FALSE,
    thealpha = 0.05,
    maxtest = 10,
    trace = FALSE
  )

  # Closed testing approach
  results_closed <- find_blocks(
    idat = idat,
    bdat = bdat,
    blockid = "blockF",
    splitfn = splitCluster,
    pfn = pOneway,
    fmla = Y1 ~ trtF | blockF,
    splitby = "hwt",
    parallel = "no",
    use_closed_testing = TRUE,
    closed_testing_method = "simes",
    thealpha = 0.05,
    maxtest = 10,
    trace = FALSE
  )

  # Both should produce valid results
  expect_true(is.list(results_traditional))
  expect_true(is.list(results_closed))

  # Check that both have same structure for comparison
  expect_true("bdat" %in% names(results_traditional))
  expect_true("bdat" %in% names(results_closed))
  expect_true("node_dat" %in% names(results_traditional))
  expect_true("node_dat" %in% names(results_closed))

  # Compare detection rates
  traditional_detections <- report_detections(results_traditional$bdat, fwer = TRUE)
  traditional_count <- sum(traditional_detections$hit, na.rm = TRUE)

  expect_true(is.numeric(traditional_count))
  expect_true(traditional_count >= 0)

  # If closed testing was applied, check its results
  if ("closed_testing_reject" %in% names(results_closed$node_dat)) {
    closed_count <- sum(results_closed$node_dat$closed_testing_reject, na.rm = TRUE)
    expect_true(is.numeric(closed_count))
    expect_true(closed_count >= 0)
  }
})

test_that("find_blocks parameters validation for advanced methods", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  # Test that valid parameters work
  expect_no_error({
    find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pOneway,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_closed_testing = FALSE,
      use_evalues = FALSE,
      closed_testing_method = "simes",
      evalue_wealth_rule = "kelly",
      maxtest = 5,
      trace = FALSE
    )
  })

  # Test boolean parameters
  expect_no_error({
    find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pOneway,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_closed_testing = TRUE,
      use_evalues = FALSE,
      maxtest = 5,
      trace = FALSE
    )
  })
})

test_that("advanced methods work with different splitting functions", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

  data(example_dat, package = "manytestsr")

  idat <- data.table::as.data.table(example_dat)
  bdat <- idat %>%
    dplyr::group_by(blockF) %>%
    dplyr::summarize(
      nb = dplyr::n(),
      pb = mean(trt),
      hwt = (nb / nrow(idat)) * (pb * (1 - pb)),
      place = unique(place),
      year = unique(year),
      place_year_block = factor(unique(place_year_block)),
      .groups = "drop"
    ) %>%
    data.table::as.data.table()

  splitting_functions <- list(
    list(fn = splitCluster, splitby = "hwt"),
    list(fn = splitEqualApprox, splitby = "hwt"),
    list(fn = splitSpecifiedFactor, splitby = "place_year_block")
  )

  for (i in seq_along(splitting_functions)) {
    split_info <- splitting_functions[[i]]

    expect_no_error({
      result <- find_blocks(
        idat = idat,
        bdat = bdat,
        blockid = "blockF",
        splitfn = split_info$fn,
        pfn = pIndepDist,
        fmla = Y1 ~ trtF | blockF,
        splitby = split_info$splitby,
        parallel = "no",
        use_closed_testing = TRUE,
        closed_testing_method = "simes",
        maxtest = 8,
        trace = FALSE
      )
    })

    expect_true(is.list(result))
    expect_true("node_dat" %in% names(result))
  }
})

test_that("advanced methods work with different test functions", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  test_functions <- list(pOneway, pIndepDist, pWilcox)

  for (pfn in test_functions) {
    expect_no_error({
      result <- find_blocks(
        idat = idat,
        bdat = bdat,
        blockid = "blockF",
        splitfn = splitCluster,
        pfn = pfn,
        fmla = Y1 ~ trtF | blockF,
        splitby = "hwt",
        parallel = "no",
        use_closed_testing = TRUE,
        closed_testing_method = "simes",
        maxtest = 6,
        trace = FALSE
      )
    })

    expect_true(is.list(result))
    expect_true("node_dat" %in% names(result))
  }
})

test_that("sequential FDR methods work with advanced approaches", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  # Test closed testing with alpha investing
  expect_no_error({
    result_alpha_investing <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pIndepDist,
      alphafn = alpha_investing,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_closed_testing = TRUE,
      closed_testing_method = "simes",
      thealpha = 0.05,
      thew0 = 0.049,
      maxtest = 8,
      trace = FALSE
    )
  })

  # Test with SAFFRON
  expect_no_error({
    result_saffron <- find_blocks(
      idat = idat,
      bdat = bdat,
      blockid = "blockF",
      splitfn = splitCluster,
      pfn = pIndepDist,
      alphafn = alpha_saffron,
      fmla = Y1 ~ trtF | blockF,
      splitby = "hwt",
      parallel = "no",
      use_closed_testing = TRUE,
      closed_testing_method = "simes",
      thealpha = 0.05,
      thew0 = 0.049,
      maxtest = 8,
      trace = FALSE
    )
  })

  expect_true(is.list(result_alpha_investing))
  expect_true(is.list(result_saffron))
})

test_that("performance with advanced methods is reasonable", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("dplyr")

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

  # Test that advanced methods complete in reasonable time
  start_time <- Sys.time()

  result <- find_blocks(
    idat = idat,
    bdat = bdat,
    blockid = "blockF",
    splitfn = splitCluster,
    pfn = pOneway,
    fmla = Y1 ~ trtF | blockF,
    splitby = "hwt",
    parallel = "no",
    use_closed_testing = TRUE,
    closed_testing_method = "simes",
    maxtest = 15, # Reasonable limit for testing
    trace = FALSE
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Should complete within 30 seconds on most systems
  expect_true(elapsed < 30,
    info = paste("Advanced methods took", round(elapsed, 2), "seconds")
  )

  expect_true(is.list(result))
})
