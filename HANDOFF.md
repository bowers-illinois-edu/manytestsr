# Handoff Summary

**Date:** 2026-02-24 **Branch:** `main` **Last commit:** `f5ec1df` —
Update README and fix LICENSE.md to match MIT declaration

------------------------------------------------------------------------

## 1. Key Decisions Made

- **License is MIT.** The `DESCRIPTION` and `LICENSE` files already
  declared MIT. The `LICENSE.md` file contained GPL-2 text (likely from
  a project template). Replaced it with the MIT license text to resolve
  the inconsistency. User confirmed MIT is the intended license.
- **README rewritten for public-facing clarity.** Added a
  development-stage warning, installation instructions, a public TODO
  checklist with done/open items, a key functions table, and updated
  development commands to reflect the Makefile workflow.

## 2. Files Changed and Why

### `README.md`

Full rewrite. The previous README had a brief warning, a three-item TODO
list, outdated `devtools::` commands, and an implementation note with a
typo. Replaced with:

- Package description with citation link to Bowers and Chen (2020)
- Bold development-stage warning placed before installation instructions
- [`remotes::install_github()`](https://remotes.r-lib.org/reference/install_github.html)
  installation instructions with C++ compiler note
- Public TODO checklist organized by category (alpha adjustment, local
  p-value adjustment, testing, documentation) using GitHub checkbox
  syntax
- Key functions reference table
- `make` commands for development
- Cleaned-up implementation notes on the C++ code paths

### `LICENSE.md`

Replaced GPL-2 full text with MIT license text (Copyright 2020 Jake
Bowers). Now consistent with `DESCRIPTION`
(`License: MIT + file LICENSE`) and the `LICENSE` file (YEAR/COPYRIGHT
HOLDER format).

## 3. Current Blockers or Open Questions

- **None.** All changes are committed and pushed to `origin/main`.

## 4. Important Context to Preserve

- **`.Rbuildignore` is regex-based.** Every line is a Perl regular
  expression, not a glob or literal path. Unanchored patterns can
  silently exclude files deep in the directory tree. When adding
  entries, always anchor with `^` and escape dots with `\.`.
- **License files:** Three files must stay in sync — `DESCRIPTION`
  (declares `MIT + file LICENSE`), `LICENSE` (short YEAR/COPYRIGHT
  HOLDER format for CRAN), and `LICENSE.md` (full text for GitHub
  display).
- **Development workflow:** `make document`, `make test`, `make check`
  are the core commands. The project uses `renv` for dependency
  management. Tests exclude profiling by default. See `CLAUDE.md` for
  full conventions.
- **Package version:** Currently `0.0.4.1002`. Recent additions include
  tree-based adaptive alpha (`alpha_adaptive_tree`), branch-pruning
  alpha (`alpha_adaptive_tree_pruned`), and error-load diagnostics
  (`compute_error_load`). See `NEWS.md` for full changelog.

## 5. What’s Done vs. What Remains

### Done

- README rewritten with development warning, TODO list, install
  instructions, and key functions
- LICENSE.md corrected from GPL-2 to MIT
- All changes committed and pushed to `origin/main`

### Remains

- Nothing from this session is outstanding
- The public TODO items in README.md track ongoing package development
  work (not session tasks)
