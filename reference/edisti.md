# Outcome e-distances between treatment arms

Outcome e-distances between treatment arms

## Usage

``` r
edisti(x, Z)
```

## Arguments

- x:

  Is a numeric vector (the outcome variable)

- Z:

  is a binary numeric vector or a factor vector with only two values

## Value

Vector with the individual level components of the energy distances
between each unit and the units in the control condition

## Examples

``` r
# Example with treatment and control groups
outcome <- c(2.1, 5.3, 1.8, 7.2, 3.4, 4.6, 6.1, 2.8)
treatment <- c(0, 1, 0, 1, 0, 1, 1, 0) # Binary treatment indicator

# Compute energy distances
e_dists <- edisti(outcome, treatment)
print(e_dists)
#>      1      2      3      4      5      6      7      8 
#> 1.5625 0.9625 1.6375 1.6375 0.7625 0.4375 1.3625 1.2125 

if (FALSE) { # \dontrun{
# With factor treatment variable
treatment_factor <- factor(treatment, labels = c("Control", "Treatment"))
e_dists2 <- edisti(outcome, treatment_factor)
print(e_dists2)
} # }
```
