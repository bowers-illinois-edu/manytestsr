# Return detected blocks plus info

Given the results of the splitting and testing algorithm, report on the
blocks where the null of no effects could be rejected at level alpha.
Currently calculates rejections using an FWER style criteria (p of a
node = max of all previous nodes) if the final alphas are all the same
as the scalar alpha OR if fwer=TRUE.

## Usage

``` r
report_detections(
  orig_res,
  fwer = TRUE,
  alpha = 0.05,
  only_hits = FALSE,
  blockid = "blockF"
)
```

## Arguments

- orig_res:

  results data.table output from the
  [`find_blocks`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)
  function.

- fwer:

  (default is TRUE) means that a block is detected (or not) using the
  maximum p-value associated with the block (or the groups containing
  that block). fwer=FALSE to detect blocks (or groups of blocks) using
  FDR control.

- alpha:

  Is the false positive rate used for detecting an effect if it is
  constant (i.e. not an FDR-style approach).

- only_hits:

  (default FALSE) returns only the detected blocks instead of all of
  them

- blockid:

  Name of block variable (the blocking variable is a factor)

## Value

A data.table adding columns to the `res` data.table: `hit` (TRUE when
the block is detected either by its own test or by an unattributed group
rejection; never NA); `hit_type` ("single" = the block's own singleton
test rejected; "group" = the block's family parent rejected while no
child of that parent rejected, so the effect is localized to the parent
but not attributed to specific blocks; "none" otherwise – including
blocks under a rejected parent that a sibling's individual rejection
already explains); `group_p` (for "group" blocks, the rejecting parent's
p-value; NA otherwise); and `fin_depth` (the block's own final tested
depth). `pfinalb` from
[`find_blocks`](https://bowers-illinois-edu.github.io/manytestsr/reference/find_blocks.md)
is the running maximum of p-values along the block's tested path, so for
"group" blocks it shows the non-significant child-level p while
`group_p` shows the parent rejection that constitutes the finding.

## Examples

``` r
if (FALSE) { # \dontrun{
# Use example data and run find_blocks
data(example_dat, package = "manytestsr")
library(data.table)
library(dplyr)

# Create block-level dataset
example_bdat <- example_dat %>%
  group_by(blockF) %>%
  summarize(
    nb = n(),
    pb = mean(trt),
    hwt = (nb / nrow(example_dat)) * (pb * (1 - pb)),
    .groups = "drop"
  ) %>%
  as.data.table()

# Run find_blocks
results <- find_blocks(
  idat = example_dat,
  bdat = example_bdat,
  blockid = "blockF",
  splitfn = splitCluster,
  pfn = pOneway,
  fmla = Y1 ~ trtF | blockF,
  parallel = "no"
)

# Report detections using FWER control
detections_fwer <- report_detections(results$bdat, fwer = TRUE, alpha = 0.05)
head(detections_fwer[, .(blockF, hit, pfinalb)])

# Report only significant blocks
hits_only <- report_detections(results$bdat, fwer = TRUE, only_hits = TRUE)
print(hits_only)
} # }
```
