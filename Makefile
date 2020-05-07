# load.R fixes a bug with devtool's `help` to enable `help` on
# functions in this package, as well as loading the package
LOAD=R_PROFILE=load.R
RCMD=R_LIBS=.local R -q 

interactive:
	@$(LOAD) $(RCMD) --no-save

.devtools:
	@$(RCMD) -e "devtools:::$(FUNC)($(DEVTOOLSARG))"

DEVTOOLSARG=
dependencies: FUNC=install_deps
dependencies: DEVTOOLSARG=dependencies=TRUE
test: FUNC=test
check: FUNC=check
document: FUNC=document
# vignette: FUNC=build_vignettes # To be renabled if we add vignettes
# clean-vignette: FUNC=clean_vignettes
build: FUNC=build
dependencies test check document build: .devtools
#vignette clean-vignette: .devtools

clean: #clean-vignette
	git clean -Xfd
