## Misc tests of different helper functions



#### Testing the C++ functions

test_that("avg_rank_arma matches base::rank for random data", {
  set.seed(42)
  x <- rnorm(100)
  expect_equal(avg_rank_arma(x),
    rank(x, ties.method = "average"),
    tolerance = 1e-12
  )
})

test_that("avg_rank_arma handles ties correctly", {
  x <- c(3, 3, 7, 7, 7, 10)
  expect_equal(avg_rank_arma(x),
    rank(x, ties.method = "average"),
    tolerance = 1e-12
  )
})

test_that("avg_rank_arma works for length-0 and length-1 vectors", {
  expect_equal(avg_rank_arma(numeric(0)), numeric(0))
  expect_equal(avg_rank_arma(42), 1)
})


# mean absolute distance from each unit i to every other unit j
mean_dist_R <- function(x) {
  sapply(seq_along(x), function(i) mean(abs(x[i] - x[-i])))
}

# max absolute distance for each unit i
max_dist_R <- function(x) {
  sapply(seq_along(x), function(i) max(abs(x[i] - x[-i])))
}


test_that("mean and max distances match R reference", {
  set.seed(20250619)
  x <- rnorm(20) # no ties
  res <- fast_dists_and_trans_new(x)

  expect_equal(res$mean_dist, mean_dist_R(x), tolerance = 1e-12)
  expect_equal(res$max_dist, max_dist_R(x), tolerance = 1e-12)

  rnk <- rank(x, ties.method = "average")
  expect_equal(res$mean_rank_dist, mean_dist_R(rnk), tolerance = 1e-12)
})

test_that("ties are handled exactly like rank(..., ties = 'average')", {
  x <- c(3, 3, 7, 7, 7, 10)
  res <- fast_dists_and_trans_new(x)

  rank_ref <- rank(x, ties.method = "average")
  expect_equal(res$rank, rank_ref, tolerance = 1e-12)

  expect_equal(res$mean_dist, mean_dist_R(x), tolerance = 1e-12)
  expect_equal(res$max_dist, max_dist_R(x), tolerance = 1e-12)
})

test_that("works for minimal length n = 2", {
  x <- c(5, 11)
  res <- fast_dists_and_trans_new(x)

  expect_equal(res$mean_dist, c(6, 6))
  expect_equal(res$mean_rank_dist, c(1, 1)) # ranks are 1 and 2 → |1-2| = 1
  expect_equal(res$max_dist, c(6, 6))
})

test_that("tanh transform is correct", {
  x <- seq(-2, 2, length.out = 5)
  res <- fast_dists_and_trans_new(x)
  expect_equal(res$tanh, tanh(x), tolerance = 1e-12)
})

test_that("OMP and serial results are identical", {
  set.seed(1)
  x <- rnorm(101)
  expect_equal(fast_dists_and_trans_new_omp(x, threads = 2),
    fast_dists_and_trans_new(x),
    tolerance = 1e-12
  )
})
