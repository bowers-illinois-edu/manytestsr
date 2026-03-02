## Tests for pPolyRank() — polynomial rank score test via coin

test_that("pPolyRank returns a numeric scalar between 0 and 1", {
  p <- pPolyRank(as.data.table(idat), Y ~ ZF | bF)
  expect_type(p, "double")
  expect_length(p, 1)
  expect_true(p >= 0 && p <= 1)
})

test_that("pPolyRank rejects on data with known effects", {
  # Yhomog has strong uniform positive effects
  p <- pPolyRank(as.data.table(idat), Yhomog ~ ZF | bF)
  expect_lt(p, 0.001)
})

test_that("pPolyRank does not reject on null data", {
  p_null <- pPolyRank(as.data.table(idat), Ynull ~ ZF | bF)
  expect_gt(p_null, 0.05)
})

test_that("pPolyRank returns 1 for constant outcome", {
  idat_const <- copy(as.data.table(idat))
  idat_const[, Yconst := 5]
  p <- pPolyRank(idat_const, Yconst ~ ZF | bF)
  expect_equal(p, 1)
})

test_that("pPolyRank works without block variable", {
  p <- pPolyRank(as.data.table(idat), Yhomog ~ ZF)
  expect_type(p, "double")
  expect_lt(p, 0.001)
})

test_that("pPolyRank works with teststat = maximum", {
  p <- pPolyRank(as.data.table(idat), Yhomog ~ ZF | bF, teststat = "maximum")
  expect_type(p, "double")
  expect_lt(p, 0.001)
})

test_that("pPolyRank works with a single r value", {
  # With r=2 this should behave like a Wilcoxon-type test
  p <- pPolyRank(as.data.table(idat), Yhomog ~ ZF | bF, r_vec = 2)
  expect_type(p, "double")
  expect_lt(p, 0.001)
})
