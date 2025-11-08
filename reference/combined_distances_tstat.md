# The combined distances test statistic

The combined distances test statistic

## Usage

``` r
combined_distances_tstat(
  fmla = Y ~ trtF | blockF,
  distfn = fast_dists_and_trans_hybrid,
  dat
)
```

## Arguments

- fmla:

  A formula specifying the response, treatment, and optionally block
  variables (e.g., Y ~ trtF \| blockF)

- distfn:

  A distance function to compute distances (default:
  fast_dists_and_trans_hybrid)

- dat:

  A data.table containing the variables specified in the formula

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires dat object to be defined in environment
# tstat_dt <- combined_distances_tstat()

splitCluster("blockF", x)
## clus <- tryCatch(KMeans_rcpp(as.matrix(x), clusters = 2, num_init = 2,
##        initializer = "optimal_init")$clusters, error = function(e) {
##    kmeans(x, centers = 2)$cluster})

## Approach 2:
# clus <- Ckmeans.1d.dp(x, k = 2)$cluster
# group <- factor(as.numeric(clus == 1))
#
} # }
```
