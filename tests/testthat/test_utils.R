## Tests for the internal which_are_constant() helper
## This function replaces dataPreparation::which_are_constant()

test_that("which_are_constant finds constant columns", {
  df <- data.frame(
    a = c(1, 1, 1),
    b = c(1, 2, 3),
    c = c(5, 5, 5)
  )
  result <- which_are_constant(df)
  expect_true(1 %in% result) # column 'a' is constant
  expect_true(3 %in% result) # column 'c' is constant
  expect_false(2 %in% result) # column 'b' varies
  expect_equal(names(result), c("a", "c"))
})

test_that("which_are_constant returns integer(0) when no columns are constant", {
  df <- data.frame(
    a = c(1, 2, 3),
    b = c(4, 5, 6)
  )
  result <- which_are_constant(df)
  expect_length(result, 0)
  expect_type(result, "integer")
})

test_that("which_are_constant returns all indices when all columns are constant", {
  df <- data.frame(
    a = c(7, 7, 7),
    b = c(3, 3, 3)
  )
  result <- which_are_constant(df)
  expect_equal(length(result), 2)
  expect_equal(as.integer(result), c(1L, 2L))
})

test_that("which_are_constant handles single-row data frames", {
  df <- data.frame(a = 1, b = 2)
  result <- which_are_constant(df)
  # A single-row data frame has all constant columns (fewer than 2 unique values)
  expect_equal(length(result), 2)
})

test_that("which_are_constant handles NA values", {
  # A column with only NAs has fewer than 2 unique values
  df <- data.frame(
    a = c(NA, NA, NA),
    b = c(1, NA, 3)
  )
  result <- which_are_constant(df)
  expect_true(1 %in% result) # column 'a' is all NA
  expect_false(2 %in% result) # column 'b' has multiple unique values (including NA)
})

test_that("which_are_constant matches dataPreparation behavior on test data", {
  skip_if_not_installed("dataPreparation")
  # Verify our internal helper matches the original on realistic data
  df <- data.frame(
    outcome = c(1.2, 3.4, 5.6, 7.8),
    constant_col = c(0, 0, 0, 0),
    rank_col = c(1, 2, 3, 4),
    another_constant = c(NA, NA, NA, NA)
  )
  our_result <- which_are_constant(df)
  orig_result <- dataPreparation::which_are_constant(df, verbose = FALSE)
  # dataPreparation returns unnamed indices; ours are named. Compare values only.
  expect_equal(unname(our_result), unname(orig_result))
})
