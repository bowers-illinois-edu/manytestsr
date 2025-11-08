# Alpha adjustment function: Alpha Investing

These function accept a vector of p-values and vector indicating batch
of p-values It produces a vector of alpha values.

## Usage

``` r
alpha_investing(pval, batch, nodesize, thealpha = 0.05, thew0 = 0.05 - 0.001)
```

## Arguments

- pval:

  Is a numeric vector of p-values

- batch:

  A character or factor or even numeric vector indicating grouping among
  the p-values and the order of the p-values. The first element pval and
  first element of batch should be the first p-value produced in the
  tree and the subsequent p-values should be produced on the children of
  that root (or children of the children of the root, etc.)

- nodesize:

  Is vector indicating the information used to create each p-value. For
  example, it could be the number of observations, or a weighted
  version.

- thealpha:

  Is the overall error rate for a given test

- thew0:

  Is the starting "wealth" of the alpha investing procedure

## Value

A vector of alpha values

## Examples

``` r
# Example with hierarchical testing scenario
# Simulate p-values from different levels of a tree
pvals <- c(0.01, 0.04, 0.12, 0.08, 0.15, 0.02)
batches <- c(1, 2, 2, 3, 3, 3) # Tree depth levels
node_sizes <- c(100, 50, 50, 25, 25, 25) # Sample sizes at each node

# Apply alpha investing procedure
alpha_vals <- alpha_investing(pvals, batches, node_sizes, thealpha = 0.05)

# Compare p-values to adjusted alpha levels
data.frame(
  pval = pvals,
  batch = batches,
  alpha = alpha_vals,
  significant = pvals <= alpha_vals
)
#>   pval batch       alpha significant
#> 1 0.01     1 0.020987117        TRUE
#> 2 0.04     2 0.021406257       FALSE
#> 3 0.12     2 0.007164201       FALSE
#> 4 0.08     3 0.003757589       FALSE
#> 5 0.15     3 0.002374706       FALSE
#> 6 0.02     3 0.001662890       FALSE
```
