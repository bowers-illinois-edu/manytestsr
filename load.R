## Use devtools to load most recent version
library(utils)
library("devtools")
devtools:::load_all(export_all=FALSE)
setDTthreads(1)

## Dealing with the "no visiblle binding for global variable" problem
##utils::globalVariables(c("prob", "section", "y"))

## `utils` is loaded first to ensure proper ordering in the search
## stack, so that devtool's versions of `help` and `?` mask `utils`.
