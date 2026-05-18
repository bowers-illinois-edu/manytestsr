# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package.
#
# First-time setup: `make dependencies` (installs CRAN deps, GitHub Remotes,
# and LinkingTo deps via pak). After that, `make check` should run cleanly.
LOAD=R_PROFILE=load.R
RCMD=R -q

.PHONY:interactive
interactive:
	@$(LOAD) $(RCMD) --no-save

.PHONY:.devtools
.devtools:
	@$(RCMD) -e "devtools:::$(FUNC)($(DEVTOOLSARG))"

DEVTOOLSARG=

# Install all declared dependencies (Depends, Imports, LinkingTo, Suggests)
# plus GitHub Remotes via pak. pak handles Remotes natively and resolves
# system requirements better than devtools::install_deps.
.PHONY: dependencies
dependencies:
	@$(RCMD) -e "if (!requireNamespace('pak', quietly = TRUE)) install.packages('pak', repos = 'https://cloud.r-project.org'); pak::pkg_install('.', dependencies = TRUE, upgrade = FALSE, ask = FALSE)"

.PHONY: test
test: FUNC=test

# Standard package check. Vignettes are built (requires Quarto CLI).
# To skip vignettes, use `make check-fast`.
.PHONY: check
check: FUNC=check

.PHONY: check-fast
check-fast: FUNC=check
check-fast: DEVTOOLSARG=args = c('--no-build-vignettes', '--no-manual'), build_args = '--no-build-vignettes', vignettes = FALSE, manual = FALSE

.PHONY: check-cran
check-cran: FUNC=check
check-cran: DEVTOOLSARG=cran = TRUE, remote = TRUE, manual = TRUE

.PHONY: document
document: FUNC=document
# vignette: FUNC=build_vignettes # To be renabled if we add vignettes
# clean-vignette: FUNC=clean_vignettes
#
.PHONY: build
build: FUNC=build
test check check-fast check-cran document build: .devtools
#vignette clean-vignette: .devtools

.PHONY:clean
clean: #clean-vignette
	git clean -Xfd

.PHONY:spell-check-DESCRIPTION
spell-check-DESCRIPTION:
	aspell -c DESCRIPTION --personal=NULL
