# Example Data for Block-Randomized Experiment

A simulated dataset demonstrating a block-randomized experiment with two
outcome variables. This dataset is used throughout the package examples
to illustrate the use of various testing and splitting procedures.

## Usage

``` r
example_dat
```

## Format

A data.table with 1268 rows and 9 variables:

- id:

  Integer. Unique identifier for each unit.

- year:

  Integer. Year indicator (1 or 3).

- trt:

  Integer. Binary treatment indicator (0 = control, 1 = treated).

- Y1:

  Numeric. First outcome variable.

- Y2:

  Numeric. Second outcome variable.

- trtF:

  Factor. Treatment indicator as a factor with levels "0" and "1".

- place_year_block:

  Character. Combined identifier for place, year, and block.

- place:

  Character. Place/location identifier (A, B, C, etc.).

- blockF:

  Factor. Block identifier with 44 levels (B080 through B123).

## Source

Simulated data for demonstration purposes

## Examples

``` r
data(example_dat, package = "manytestsr")
head(example_dat)
#>       id  year   trt    Y1    Y2   trtF place_year_block  place blockF
#>    <int> <int> <int> <num> <num> <fctr>           <char> <char> <fctr>
#> 1:     1     1     0     0     0      0         A.1.B082      A   B082
#> 2:     2     3     0     0    12      0         B.3.B094      B   B094
#> 3:     3     1     0     0     0      0         C.1.B097      C   B097
#> 4:     4     1     0     6     0      0         C.1.B097      C   B097
#> 5:     5     1     0     7    11      0         B.1.B089      B   B089
#> 6:     6     1     1     0     0      1         A.1.B080      A   B080
table(example_dat$trt)
#> 
#>   0   1 
#> 439 829 
table(example_dat$blockF)
#> 
#> B080 B081 B082 B083 B084 B085 B086 B087 B088 B089 B090 B091 B092 B093 B094 B095 
#>  129   68   56    8    9   53  167   60   23   62    4    6    2    5   42   14 
#> B096 B097 B098 B099 B100 B101 B102 B103 B104 B105 B106 B107 B108 B109 B110 B111 
#>   14   52   20   12    9    4   41   23   21    3   34    3    3    2    4   38 
#> B112 B113 B114 B115 B116 B117 B118 B119 B120 B121 B122 B123 
#>   10    9    2   54   24   27   14   13   34   21   50   19 
```
