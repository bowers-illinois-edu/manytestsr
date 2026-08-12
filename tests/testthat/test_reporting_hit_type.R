# Tests for honest detection reporting (referee-correction Fix 3, see
# FIX_PLAN.md). Written 2026-08-11 BEFORE the implementation fix.
#
# The statistical principle: a block whose OWN test rejected and a block
# merely COVERED by a rejected group hypothesis are different findings.
# The procedure rejected the joint null at the parent; it did not separate
# the children. Reporting both as "hit" with no distinction (and, in the
# supplement's F.3-F.4 case, with pfinalb showing a non-significant 0.134)
# is what the AOAS referee called not-like-for-like. These tests demand:
#   - hit_type in {"single", "group", "none"}: "single" = own test
#     rejected; "group" = under a rejected parent whose children ALL
#     failed individually (the procedure could not separate them);
#     "none" = everything else, including the mixed-family case where a
#     sibling's individual rejection accounts for the group finding;
#   - hit is never NA (the pre-fix code returns NA for branches that stop
#     early, which sum(hit, na.rm = TRUE) callers silently drop);
#   - group coverage detected on branches stopping BEFORE the global
#     maximum depth (pre-fix group_hit reads the global max-depth p
#     column, so shallower branches silently get FALSE or NA);
#   - pfinalb equal to the running MAXIMUM of p-values on the block's
#     tested path (the weakest evidence on the path), not the deepest p;
#   - group_p carrying the rejecting parent's p-value for covered blocks.
#
# Deterministic fixture (asymptotic test branch; seed picked by margin
# search so every precondition holds with slack; validated in-test):
#   Site A: 2 cohorts x 2 blocks (depth-4 leaves), tau = 2 -> single hits.
#   Site B: 4 blocks, constant cohort field (depth-3 leaves), tau = 0.25
#           -> site rejects (p ~ .034), every block fails individually
#           (min p ~ .085) -> pure group coverage on a shallow branch.
#   Site C: 2 blocks, tau = 0 -> site fails (p ~ .45); never tested
#           deeper; must be hit = FALSE, not NA.
#   Site D: 2 blocks, tau = (1.0, 0) -> site rejects (p ~ .0077), D1
#           rejects individually (p ~ .00036), D2 fails (p ~ .94). Mixed
#           family: D1 "single", D2 "none". Because site D's p exceeds
#           D1's own p, D1's running-max pfinalb differs from its deepest
#           p -- this is the case that makes the pmax test discriminate.

library(data.table)

make_reporting_fixture <- function() {
  set.seed(1)
  block_spec <- rbind(
    data.frame(
      site = "A", cohort = rep(c("c1", "c2"), each = 2),
      block = c("A1", "A2", "A3", "A4"), n = 60, tau = 2.0
    ),
    data.frame(
      site = "B", cohort = "x", block = c("B1", "B2", "B3", "B4"),
      n = 80, tau = 0.25
    ),
    data.frame(
      site = "C", cohort = "x", block = c("C1", "C2"), n = 60, tau = 0
    ),
    data.frame(
      site = "D", cohort = "x", block = c("D1", "D2"), n = 80,
      tau = c(1.0, 0)
    )
  )
  idat <- do.call(rbind, lapply(seq_len(nrow(block_spec)), function(i) {
    s <- block_spec[i, ]
    trt <- rep(c(0L, 1L), each = s$n / 2)
    data.frame(
      blockF = s$block,
      scb = paste(s$site, s$cohort, s$block, sep = "."),
      trtF = factor(trt),
      Y1 = rnorm(s$n) + s$tau * trt
    )
  }))
  idat <- as.data.table(idat)
  idat[, blockF := factor(blockF)]
  idat[, scb := factor(scb)]
  bdat <- idat[, .(nb = .N, scb = unique(scb)), by = blockF]
  list(idat = idat, bdat = bdat)
}

fx <- make_reporting_fixture()
res <- find_blocks(
  idat = fx$idat, bdat = fx$bdat, blockid = "blockF",
  pfn = pOneway, fmla = Y1 ~ trtF | blockF,
  splitfn = splitSpecifiedFactorMulti,
  alphafn = NULL,
  splitby = "scb",
  blocksize = "nb",
  copydts = TRUE, ncores = 1, parallel = "no", trace = FALSE
)
det <- report_detections(
  res$bdat,
  fwer = TRUE, alpha = 0.05, only_hits = FALSE, blockid = "blockF"
)
setkey(det, blockF)

p_cols <- grep("^p[0-9]+$", names(res$bdat), value = TRUE)
bd <- res$bdat[order(blockF)]
path_p <- as.matrix(bd[, ..p_cols])
rownames(path_p) <- as.character(bd$blockF)
blocks_of <- function(prefix) grep(paste0("^", prefix), rownames(path_p))

test_that("fixture behaves as designed (validation, not the claim)", {
  deepest <- apply(path_p, 1, function(z) tail(na.omit(z), 1))
  # Site A: every block's own (deepest) test rejects.
  expect_true(all(deepest[blocks_of("A")] < 0.05))
  # Site B: the site gate rejected but every block's own test failed.
  expect_true(all(deepest[blocks_of("B")] > 0.05))
  expect_true(all(path_p[blocks_of("B"), "p2"] < 0.05))
  # Site C: stopped at the site level.
  expect_true(all(path_p[blocks_of("C"), "p2"] > 0.05))
  expect_true(all(is.na(deepest[blocks_of("C")]) | deepest[blocks_of("C")] > 0.05))
  # Site D: mixed family, and the parent p exceeds D1's own p so the
  # running-max test below can discriminate from deepest-p behavior.
  expect_lt(deepest["D1"], 0.001)
  expect_gt(deepest["D2"], 0.10)
  expect_gt(path_p["D1", "p2"], deepest["D1"])
  expect_lt(path_p["D1", "p2"], 0.05)
})

test_that("hit_type separates single, group-covered, and unattributed", {
  expect_true("hit_type" %in% names(det))
  ht <- setNames(as.character(det$hit_type), as.character(det$blockF))
  expect_setequal(unname(ht[paste0("A", 1:4)]), rep("single", 4))
  expect_setequal(unname(ht[paste0("B", 1:4)]), rep("group", 4))
  expect_setequal(unname(ht[c("C1", "C2")]), rep("none", 2))
  expect_equal(unname(ht["D1"]), "single")
  # D2 sits under a rejected parent, but its sibling was individually
  # separated, so the group finding is attributed to D1, not spread over
  # the family.
  expect_equal(unname(ht["D2"]), "none")
})

test_that("hit is never NA, and shallow-branch group coverage is found", {
  # Pre-fix behavior on this fixture (run-verified): site C blocks get
  # hit = NA because their branch stopped before the global max depth,
  # and site B blocks get FALSE for the same structural reason even
  # though their parent's group null was rejected.
  expect_false(any(is.na(det$hit)))
  hits <- setNames(det$hit, as.character(det$blockF))
  expect_true(all(hits[paste0("B", 1:4)]))
  expect_false(any(hits[c("C1", "C2")]))
  expect_false(unname(hits["D2"]))
})

test_that("pfinalb is the running maximum of the tested path", {
  expected <- apply(path_p, 1, max, na.rm = TRUE)
  got <- setNames(det$pfinalb, as.character(det$blockF))
  # D1 is the discriminating case: deepest p is 4e-4 but the weakest
  # evidence on its path is the site-level ~ .0077.
  expect_equal(unname(got[names(expected)]), unname(expected),
    tolerance = 1e-12
  )
})

test_that("covered blocks carry the rejecting parent's p-value", {
  expect_true("group_p" %in% names(det))
  gp <- setNames(det$group_p, as.character(det$blockF))
  b_expected <- unname(path_p["B1", "p2"])
  expect_equal(unname(gp[paste0("B", 1:4)]), rep(b_expected, 4),
    tolerance = 1e-12
  )
  expect_true(all(gp[paste0("B", 1:4)] < 0.05))
})
