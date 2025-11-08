# A set of pre-specified splits

This function does binary splits using a factor variable with dots
separating the names of the subgroups. If there are more than two
subgroups at any level, it makes one group from the largest subgroup and
another group from the rest. If there are multiple subgroups the same
size, it chooses the first subgroup by order of the levels of the
factor.

## Usage

``` r
splitSpecifiedFactor(bid, x)
```

## Arguments

- bid:

  Block id

- x:

  Is a a factor with levels like "state.district.school". The splits
  will occur from left to right depending on whether there is existing
  variation at that level

## Examples

``` r
# Create hierarchical factor for pre-specified splitting
block_ids <- c("B1", "B2", "B3", "B4", "B5", "B6")
hierarchical_factor <- factor(c(
  "StateA.District1.School1",
  "StateA.District1.School2",
  "StateA.District2.School1",
  "StateB.District1.School1",
  "StateB.District1.School2",
  "StateB.District2.School1"
))

# First split will separate by state (StateA vs StateB)
groups <- splitSpecifiedFactor(block_ids, hierarchical_factor)

# Show the grouping
data.frame(block = block_ids, hierarchy = hierarchical_factor, group = groups)
#>   block                hierarchy group
#> 1    B1 StateA.District1.School1     1
#> 2    B2 StateA.District1.School2     1
#> 3    B3 StateA.District2.School1     1
#> 4    B4 StateB.District1.School1     0
#> 5    B5 StateB.District1.School2     0
#> 6    B6 StateB.District2.School1     0
```
