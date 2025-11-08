# Fast per-unit distance summaries (scalar outcome)

Computes, for **each observation** \\i = 1,\dots,n\\ in a numeric vector
\\x\\, scalar summaries plus simple transforms, **without ever forming
the full \\n\times n\\ distance matrix**.

## Usage

``` r
fast_dists_and_trans_hybrid(x)
```

## Arguments

- x:

  Numeric vector – the scalar outcome for all units.

## Value

A named `list` with components

- `mean_dist`:

  numeric vector, length `n`.

- `mean_rank_dist`:

  numeric vector, length `n`.

- `max_dist`:

  numeric vector, length `n`.

- `rankY`:

  average ranks (mid-ranks).

- `tanhY`:

  \\\tanh(x_i)\\ values.

## Details

- **`mean_dist`** – mean absolute distance \$\$\frac{1}{n-1}\sum\_{j\neq
  i}\|x_i-x_j\|\$\$

- **`mean_rank_dist`** – same mean on the mid-rank scale; closed-form,
  no second loop.

- **`max_dist`** – maximum absolute distance
  \\\max\\\\x_i-\min(x),\\\max(x)-x_i\\\\\\

- **`rankY`** – average (mid-) rank of `x` (`ties="average"`).

- **`tanhY`** – element-wise \\\tanh(x_i)\\ shrink transform.

**Complexity**

- \\O(n \log n)\\ time (sorting + prefix sums)

- \\O(n)\\  space (only vectors of length `n`)

For `n = 10\,000` this typically runs in ≈2 ms on an Apple M-series core
with \< 0.5 MB peak RAM, much faster and far lighter than allocating the
full distance matrix.

## Examples

``` r
set.seed(1)
x <- rnorm(8)
fast_dists_and_trans_hybrid(x)
#> $mean_dist
#> [1] 0.9813777 0.7499214 1.1052377 1.6729445 0.7499214 1.0922432 0.7950418
#> [8] 0.9384107
#> 
#> $mean_rank_dist
#> [1] 2.571429 2.285714 4.000000 4.000000 2.285714 3.142857 2.571429 3.142857
#> 
#> $max_dist
#> [1] 2.221735 1.411637 2.430909 2.430909 1.265773 2.415749 1.323058 1.573953
#> 
#> $rankY
#> [1] 3 4 1 8 5 2 6 7
#> 
#> $tanhY
#> [1] -0.5556056  0.1816063 -0.6834867  0.9209551  0.3180784 -0.6753247  0.4521735
#> [8]  0.6281319
#> 

## compare to explicit distance matrix (slow / big):
dx <- abs(outer(x, x, "-"))
mean_dist_ref <- colSums(dx) / (length(x) - 1)
stopifnot(all.equal(
  fast_dists_and_trans_hybrid(x)$mean_dist,
  mean_dist_ref
))
```
