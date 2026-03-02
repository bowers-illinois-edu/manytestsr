# Internal utility functions

# Replaces dataPreparation::which_are_constant() to avoid pulling in a heavy
# dependency for a one-liner. Returns a named integer vector of column indices
# where the column has fewer than 2 unique values (i.e., is constant).
which_are_constant <- function(df, verbose = FALSE) {
  which(vapply(df, function(x) length(unique(x)) < 2, logical(1)))
}
