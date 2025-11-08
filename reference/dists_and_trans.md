# Outcome distances and transformations

Outcome distances and transformations

## Usage

``` r
dists_and_trans(x)
```

## Arguments

- x:

  Is a numeric vector (the outcome variable)

## Value

A list for inclusion in a data.table with distances between each unit
and other units as well as some transformations of those distances

## Examples

``` r
# Example with continuous outcome data
outcome <- c(2.1, 5.3, 1.8, 7.2, 3.4, 4.6)

# Compute distance-based transformations
dist_results <- dists_and_trans(outcome)

# View components
str(dist_results)
#> List of 5
#>  $ mean_dist     : num [1:6] 2.48 2.24 2.72 3.76 1.96 1.96
#>  $ mean_rank_dist: num [1:6] 2.2 2.2 3 3 1.8 1.8
#>  $ max_dist      : num [1:6] 5.1 3.5 5.4 5.4 3.8 2.8
#>  $ rankY         : num [1:6] 2 5 1 6 3 4
#>  $ tanhY         : num [1:6] 0.97 1 0.947 1 0.998 ...
print(dist_results$mean_dist) # Mean distance to all other units
#> [1] 2.48 2.24 2.72 3.76 1.96 1.96
print(dist_results$rankY) # Ranks of original values
#> [1] 2 5 1 6 3 4
print(dist_results$tanhY) # Hyperbolic tangent transformation
#> [1] 0.9704519 0.9999502 0.9468060 0.9999989 0.9977749 0.9997979
```
