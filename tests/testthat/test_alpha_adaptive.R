# Tests for adaptive alpha adjustment (Algorithm 1 from Appendix D)
#
# The adaptive alpha procedure adjusts significance levels at each tree
# level based on estimated power decay. The key statistical properties:
#
#   1. The root always gets nominal alpha (no adjustment needed).
#   2. When estimated power is high, deeper levels get stricter alpha
#      to compensate for the multiplicity of tests that power enables.
#   3. When the total error load (sum_G) is at most 1, natural gating
#      (the fact that low power makes rejection unlikely) suffices,
#      so nominal alpha is used.
#   4. The FWER guarantee holds when power estimates are conservative
#      (i.e., delta_hat is not underestimated).
#
# Note: compute_adaptive_alphas returns a named vector (names are depth
# levels as characters). We use [[ ]] for single-element extraction
# (drops names) and unname() for slice comparisons.

# ============================================================================
# Tests for compute_adaptive_alphas()
# ============================================================================

test_that("root alpha is always nominal regardless of parameters", {
  # The root level is a single test — no multiplicity adjustment needed.
  # This must hold for any combination of design parameters.
  params <- list(
    list(k = 2, delta_hat = 0.1, N_total = 50),
    list(k = 4, delta_hat = 0.5, N_total = 1000),
    list(k = 10, delta_hat = 1.0, N_total = 10000),
    list(k = 2, delta_hat = 0.01, N_total = 10)
  )
  for (p in params) {
    alphas <- compute_adaptive_alphas(
      k = p$k, delta_hat = p$delta_hat, N_total = p$N_total, thealpha = 0.05
    )
    expect_equal(alphas[[1]], 0.05,
      info = paste("k =", p$k, "delta =", p$delta_hat, "N =", p$N_total)
    )
  }
})

test_that("level 2 alpha is more stringent when root power is high", {
  # With large N and moderate effect size, power at the root is high.
  # The procedure should tighten alpha at level 2 to account for the
  # k children that high power enables us to test.
  alphas <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000, thealpha = 0.05
  )
  expect_lt(alphas[[2]], 0.05)
})

test_that("deep levels revert to nominal alpha when power decays", {
  # With moderate N and small effect, power decays at deeper levels.
  # The min(alpha, alpha/sum_pp) formula reverts to nominal alpha
  # when sum_pp < 1 (few tests expected to reach that depth).
  alphas <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.2, N_total = 200,
    max_depth = 10, thealpha = 0.05
  )
  # At some deep level, alpha should revert to nominal
  expect_true(any(unname(alphas[5:10]) == 0.05),
    info = "Expected some deep levels to revert to nominal alpha"
  )
})

test_that("larger k produces more stringent adjustment at same depth", {
  # More children per node means more tests at each level.
  # The adjustment must be stricter to control the FWER across
  # the larger number of tests.
  alphas_k2 <- compute_adaptive_alphas(
    k = 2, delta_hat = 0.5, N_total = 1000, thealpha = 0.05
  )
  alphas_k4 <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000, thealpha = 0.05
  )
  # At level 2, k=4 should be at least as stringent as k=2
  # (both have high power, so both should adjust)
  expect_lte(alphas_k4[[2]], alphas_k2[[2]])
})

test_that("larger delta_hat produces more stringent adjustment", {
  # Higher estimated power means the procedure expects more tests to
  # proceed to deeper levels, so it must tighten alpha more aggressively.
  alphas_small_d <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.2, N_total = 1000, thealpha = 0.05
  )
  alphas_large_d <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.8, N_total = 1000, thealpha = 0.05
  )
  # At level 2, larger delta should give smaller (more stringent) alpha
  expect_lte(alphas_large_d[[2]], alphas_small_d[[2]])
})

test_that("high error load means adjustment at deeper levels", {
  # When error load > 1 (high power, wide tree), the procedure must
  # adjust alphas at deeper levels to control FWER.
  alphas <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000,
    max_depth = 5, thealpha = 0.05
  )
  el <- attr(alphas, "error_load")
  expect_true(el$needs_adjustment)
  # All levels beyond root should be adjusted (less than or equal to nominal)
  expect_true(all(unname(alphas[2:5]) <= 0.05))
})

test_that("low error load means nominal alpha everywhere", {
  # When total error load <= 1 (low power), natural gating suffices
  # and every level gets nominal alpha.
  alphas <- compute_adaptive_alphas(
    k = 3, delta_hat = 0.01, N_total = 50,
    max_depth = 5, thealpha = 0.05
  )
  el <- attr(alphas, "error_load")
  expect_true(!el$needs_adjustment,
    info = paste("Expected sum_G <= 1 but got", el$sum_G))
  expect_true(all(unname(alphas) == 0.05))
})

test_that("very large N approaches Bonferroni-like correction", {
  # When N is enormous, power at every level approaches 1.
  # Then alpha_ell ≈ alpha / k^(ell-1), which is the Bonferroni
  # correction for the number of tests at that level.
  alphas <- compute_adaptive_alphas(
    k = 2, delta_hat = 0.5, N_total = 1e8,
    max_depth = 5, thealpha = 0.05
  )
  # With power ≈ 1, the formula gives alpha / (k^(ell-1) * 1) = alpha / k^(ell-1)
  for (ell in 2:5) {
    bonferroni <- 0.05 / (2^(ell - 1))
    expect_equal(alphas[[ell]], bonferroni, tolerance = 0.001,
      info = paste("level", ell)
    )
  }
})

test_that("very small N gives nominal alpha at all levels", {
  # When N is tiny, power is near zero everywhere.
  # Error load is well below 1, so natural gating suffices
  # and all levels get nominal alpha.
  alphas <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.1, N_total = 10,
    max_depth = 5, thealpha = 0.05
  )
  # Even level 2 should have nominal alpha because power is so low
  expect_equal(alphas[[2]], 0.05)
})

test_that("adjusted alpha never exceeds nominal alpha", {
  # The min(thealpha, ...) guard ensures we never give a level MORE
  # alpha than it would get without adjustment. This is important
  # for FWER control.
  alphas <- compute_adaptive_alphas(
    k = 2, delta_hat = 0.5, N_total = 1000,
    max_depth = 10, thealpha = 0.05
  )
  expect_true(all(unname(alphas) <= 0.05))
})

test_that("alpha values are always positive", {
  alphas <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000, thealpha = 0.05
  )
  expect_true(all(alphas > 0))
})

test_that("vector k is supported and differs from scalar k", {
  # When branching factor varies by level, the alpha schedule should
  # differ from a constant-k schedule.
  alphas_scalar <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000,
    max_depth = 4, thealpha = 0.05
  )
  alphas_vector <- compute_adaptive_alphas(
    k = c(4, 2, 3, 2), delta_hat = 0.5, N_total = 1000,
    max_depth = 4, thealpha = 0.05
  )
  # Root alpha is the same (always nominal)
  expect_equal(alphas_vector[[1]], alphas_scalar[[1]])
  # But deeper levels should differ because the branching factors differ
  # (level 2 has same k=4 in both, but level 3 uses k=2 vs k=4)
  expect_false(isTRUE(all.equal(alphas_vector[[3]], alphas_scalar[[3]])),
    info = "Vector and scalar k should produce different alphas at levels with different branching"
  )
})

test_that("compute_adaptive_alphas validates inputs", {
  expect_error(compute_adaptive_alphas(k = 4, delta_hat = -1, N_total = 100))
  expect_error(compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = -10))
  expect_error(compute_adaptive_alphas(k = 1, delta_hat = 0.5, N_total = 100))
  expect_error(compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 100, thealpha = 0))
  expect_error(compute_adaptive_alphas(k = 4, delta_hat = 0.5, N_total = 100, thealpha = 1))
})

test_that("known-value check for a specific parameter set", {
  # Hand-computation for k=2, N=1000, delta=0.5, alpha=0.05:
  #   z_crit = qnorm(0.975) ≈ 1.96
  #
  #   Level 1: n=1000, theta = pnorm(0.5*sqrt(1000) - 1.96)
  #            = pnorm(15.81 - 1.96) = pnorm(13.85) ≈ 1.0
  #            alpha_1 = 0.05
  #
  #   Level 2: n=500, theta = pnorm(0.5*sqrt(500) - 1.96) ≈ 1.0
  #            num_tests = 2, path_power ≈ 1.0
  #            alpha_2 = 0.05 / (2 * 1.0) = 0.025
  #
  #   Level 3: n=250, theta ≈ 1.0
  #            num_tests = 4, path_power ≈ 1.0
  #            alpha_3 = 0.05 / (4 * 1.0) = 0.0125
  alphas <- compute_adaptive_alphas(
    k = 2, delta_hat = 0.5, N_total = 1000,
    max_depth = 3, thealpha = 0.05
  )
  expect_equal(alphas[[1]], 0.05)
  expect_equal(alphas[[2]], 0.025, tolerance = 0.001)
  expect_equal(alphas[[3]], 0.0125, tolerance = 0.001)
})

# ============================================================================
# Tests for alpha_adaptive() factory
# ============================================================================

test_that("alpha_adaptive returns a function with the correct interface", {
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)
  expect_true(is.function(fn))
  # The function should accept the full alphafn interface
  args <- names(formals(fn))
  expect_true("pval" %in% args)
  expect_true("batch" %in% args)
  expect_true("nodesize" %in% args)
  expect_true("thealpha" %in% args)
  expect_true("thew0" %in% args)
  expect_true("depth" %in% args)
})

test_that("alpha_adaptive uses depth to determine alpha levels", {
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)
  # Root depth should get nominal alpha
  root_alpha <- fn(
    pval = 0.01, batch = factor("1"), nodesize = 1000,
    thealpha = 0.05, thew0 = 0.049, depth = 1L
  )
  expect_equal(unname(root_alpha), 0.05)

  # Depth 2 should get a stricter alpha (given high power at N=1000)
  depth2_alpha <- fn(
    pval = c(0.01, 0.02), batch = c("p1", "p2"),
    nodesize = c(250, 250), thealpha = 0.05, thew0 = 0.049,
    depth = c(2L, 2L)
  )
  expect_true(all(depth2_alpha < 0.05))
})

test_that("alpha_adaptive output length matches input length", {
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)
  result <- fn(
    pval = c(0.01, 0.02, 0.03),
    batch = c(1, 2, 2),
    nodesize = c(1000, 250, 250),
    thealpha = 0.05,
    thew0 = 0.049,
    depth = c(1L, 2L, 2L)
  )
  expect_length(result, 3)
})

test_that("alpha_adaptive gives consistent results with compute_adaptive_alphas", {
  # The factory's closure should produce the same alphas as calling
  # compute_adaptive_alphas directly for matching depths.
  schedule <- compute_adaptive_alphas(
    k = 4, delta_hat = 0.5, N_total = 1000,
    max_depth = 5, thealpha = 0.05
  )
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)

  for (d in 1:5) {
    result <- fn(
      pval = 0.01, batch = 1, nodesize = 100,
      thealpha = 0.05, thew0 = 0.049, depth = d
    )
    expect_equal(unname(result), schedule[[d]], info = paste("depth", d))
  }
})

test_that("alpha_adaptive handles mixed depths in a single call", {
  # find_blocks passes ancestry paths at depth >= 3, which include
  # nodes at multiple depths. The function must handle this correctly.
  fn <- alpha_adaptive(k = 2, delta_hat = 0.5, N_total = 1000)
  schedule <- compute_adaptive_alphas(
    k = 2, delta_hat = 0.5, N_total = 1000,
    max_depth = 4, thealpha = 0.05
  )

  # Simulate an ancestry path: root, depth 2, depth 3
  result <- fn(
    pval = c(0.01, 0.02, 0.03),
    batch = c(1, 2, 3),
    nodesize = c(1000, 500, 250),
    thealpha = 0.05,
    thew0 = 0.049,
    depth = c(1L, 2L, 3L)
  )
  expect_equal(unname(result[1]), schedule[[1]])
  expect_equal(unname(result[2]), schedule[[2]])
  expect_equal(unname(result[3]), schedule[[3]])
})

test_that("alpha_adaptive caches correctly and recomputes on thealpha change", {
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)

  # Two calls with same thealpha should give identical results
  res1 <- fn(pval = 0.01, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  res2 <- fn(pval = 0.99, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  expect_equal(res1, res2)

  # Different thealpha should give a different result
  res3 <- fn(pval = 0.01, batch = 1, nodesize = 250,
    thealpha = 0.10, thew0 = 0.049, depth = 2L)
  expect_false(isTRUE(all.equal(res1, res3)))
})

test_that("alpha_adaptive ignores pval (alpha depends only on tree structure)", {
  # Unlike online FDR methods, the adaptive alpha does not depend on
  # observed p-values. It depends only on design parameters and depth.
  fn <- alpha_adaptive(k = 4, delta_hat = 0.5, N_total = 1000)

  res_small_p <- fn(pval = 0.001, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  res_large_p <- fn(pval = 0.99, batch = 1, nodesize = 250,
    thealpha = 0.05, thew0 = 0.049, depth = 2L)
  expect_equal(res_small_p, res_large_p)
})

test_that("alpha_adaptive validates factory inputs", {
  expect_error(alpha_adaptive(k = 1, delta_hat = 0.5, N_total = 1000))
  expect_error(alpha_adaptive(k = 4, delta_hat = -0.5, N_total = 1000))
  expect_error(alpha_adaptive(k = 4, delta_hat = 0.5, N_total = -100))
})
