# Referee-Correction Fix Plan (2026-08-11)

Context: AOAS rejected the manytests paper (AOAS2606-025) after Referee
1 worked through this package and found implementation errors; we
verified every claim (see
`submissions/aoas/2026-06-13-initial/reviews/TRIAGE.md` and
`DECISION_MEMO.md` in `~/repos/manytests-paper`). Scope of THIS plan:
fix the package enough that it no longer misleads a user, ahead of the
paper-strategy decisions. GitHub-only distribution; the likely next
users are referees rerunning the package. Current state tagged
`pre-referee-fixes` (v0.0.4.1007) so submitted numbers stay
reproducible.

Process: tests first (this commit), Jake reviews, then implementation
(version bump to 0.0.4.1008), `make document`, `make test`,
`make check`.

## Fix 1: error-load computation no longer gives a false all-clear

All three defects bias toward understating the load, so
`needs_adjustment = FALSE` is returned when the truth is “adjust”:

- `.error_load_from_tree` (R/alpha_adaptive.R:170-246) computes
  `G[d] = sum(path_power * theta)` from depth 1. Definition 2 in the
  paper is `G_ell = sum over nodes at depth ell of path power` (proper
  ancestors only), summed from depth 2. Fix the formula; the
  `needs_adjustment` gate then automatically agrees with the
  `G_by_depth` denominators in `compute_adaptive_alphas_tree` (:544).
  Same fix in parametric `compute_error_load` (:136) and
  `compute_adaptive_alphas` (:398 already correct – only the gate at
  :364 changes source).
- One-tailed power `pnorm(delta_hat*sqrt(n) - z)` (:180, :127, :395)
  becomes two-tailed
  `pnorm(delta_hat*sqrt(n) - z) + pnorm(-delta_hat*sqrt(n) - z)`, so
  theta at delta = 0 equals the size alpha, not alpha/2.
- Units guard where nodesize enters the power formula: error if
  `any(nodesize < 1)` with a message naming the hwt trap (find_blocks’
  default `blocksize = "hwt"` produces node sizes on the unit interval;
  power at a “sample size” of 0.25 is meaningless). The guard lives in
  the error-load/power path only – `blocksize` remains free for its
  legitimate splitting/weighting roles.
- Delete the false comment at R/alpha_adaptive.R:383-388 claiming the
  two formulas agree; document that nodesize must be a headcount for all
  `alpha_adaptive*`/`compute_error_load` uses (roxygen + vignettes).

Deliberately NOT here: evaluating the load on the prespecified design
tree rather than a find_blocks node_dat (needs a tree-builder API;
paper-repo decision 6), and the ex-ante delta_hat rule. Both wait.

## Fix 2: withdraw the invalid pruned-schedule guarantee

The pruned load-denominator schedule (alpha_ell = w_ell\*alpha/D_ell)
and the switching rule implement a theorem we have now proven FALSE by
exact counterexample (FWER 0.063 at alpha = 0.05 with all hypotheses
holding; the paper’s own Job Corps sim shows 0.058). Immediate action:
[`alpha_adaptive_tree_pruned()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree_pruned.md)
issues a warning on creation stating that the strong-FWER guarantee does
not hold for the pruned-load denominator or the switching rule, and
pointing to
[`alpha_adaptive_tree()`](https://bowers-illinois-edu.github.io/manytestsr/reference/alpha_adaptive_tree.md)
(static budget weights, Theorem 4, sound). Behavior is otherwise
unchanged (warning, not removal), so old pipelines still run and the
submitted numbers remain reproducible. The replacement schedule (count
denominator or inheritance reformulation) waits on the theory-route
decision (paper-repo DECISION_MEMO decision 7).

Existing tests that deliberately pin the invalid behavior must be
rewritten at implementation time, not preserved:
`test_switching_corollary.R:223-250` (backward-compat contract),
`test_alpha_adaptive_pruned.R:191-210, :321-357` (exact-equality pins).

## Fix 3: detection reporting no longer mislabels

- [`report_detections()`](https://bowers-illinois-edu.github.io/manytestsr/reference/report_detections.md)
  gains `hit_type` in {“single”, “group”, “none”}: “single” = the
  block’s own test rejected; “group” = covered by a rejected group only.
  Pure column addition.
- `pfinalb` becomes the running maximum along the path
  (R/find_blocks.R:519; the pmax version is already there, commented
  out), so a covered block shows a p-value that is honest about what was
  and was not rejected. Audit verified traversal is unchanged (line 526
  is equivalent under constant alpha; :578 never reads pfinalb).
- Bug beyond the referee’s report (run-verified by audit): `group_hit`
  (R/reporting.R:122-127) reads the GLOBAL max-depth p column, so group
  coverage is detected only for branches reaching the global maximum
  depth; branches whose testing stopped earlier get FALSE or NA. Fix to
  use each branch’s own final depth.

## Fix 4: small things that ride along

- `coin::approximate(object, ...)` in R/pval_fns.R:137 passes a stray
  positional `object` that binds to coin’s `cl` argument (all other
  formals are matched by name); harmless under `parallel = "no"` because
  the promise is never forced, breaks under `parallel = "snow"`. Drop
  the stray argument. (See also FIXME_pval_permutation_branch.md.)
- Documentation: replace every “exact” claim about the test reference
  distribution with the true description (Strasser-Weber conditional
  moments, asymptotic shape above `simthresh`, Monte Carlo permutation
  below; no coin::exact path exists). Realized-size study waits for the
  paper rebuild.
- Lean README in the paper repo overstates formalization coverage of the
  pruning theorem; corrected there, not here (noted for completeness).

## Test plan (written first; all fail against current code by design)

- `test_error_load_definition2.R`: two-tailed size at delta = 0;
  Definition-2 value on a hand-computed fixture; root-only tree has
  sum_G = 0; gate consistent with denominators near the boundary; units
  guard on weight-scale nodesize.
- `test_pruned_schedule_guarantee.R`: warning on creation for the
  pruned-load schedule and the switching path; values still computed
  (backward compatibility).
- `test_reporting_hit_type.R`: hit_type semantics on a deterministic
  unbalanced fixture; pfinalb running-max monotonicity; group coverage
  detected on a branch that stops before the global max depth.
- `test_pval_small_n_branch.R`: small-sample branch sanity under
  `parallel = "no"` (guards the Fix-4 regression; the snow-path fix is
  asserted implicitly by the argument no longer existing).

## Environment issue found during this work (pre-existing, NOT ours)

Six test files crash R outright (“irrecoverable exception”) on this
machine under devtools::load_all + testthat, VERIFIED IDENTICAL ON THE
UNEDITED CODE for comprehensive_integration, reporting_plotting,
find_blocks_advanced, pval_fns, and misc (dist_fns crashes in the same
C++ code paths; the referee fixes touch no C++). All six exercise the
compiled distance-function layer (Rcpp/RcppArmadillo/OpenMP). Likely a
stale compilation against the current R/toolchain – try
`make clean; make build` or reinstalling with fresh compilation. Until
resolved, `make check` cannot pass on this machine for reasons unrelated
to these fixes.
