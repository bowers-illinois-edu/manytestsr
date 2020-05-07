
## This finds the individual blocks, but not the groups
# res3a <- res3long[!is.na(p), .(bF = unique(bF), maxp = max(p), hit = (max(p) <= a[p == max(p)])), by = biggrp]
# res3[, maxp := max(unlist(.SD), na.rm = TRUE), .SDcols = pnms, by = seq_len(nrow(res3))]
# res3[, maxdepth := find_maxdepth(unlist(.SD)), .SDcols = pnms, by = seq_len(nrow(res3))]


# res3lst <- split(res3, by = "biggrp")
# res3lst[1:2]
# res1lst <- split(res1, by = "biggrp")


find_maxdepth <- function(x) {
  isnax <- is.na(x)
  if (all(!isnax)) {
    maxdepth <- length(x)
  } else {
    maxdepth <- min(which(isnax)) - 1
  }
  return(maxdepth)
}


## Assessing the permutation based edist test (look into writing a custom function for this eventually)
## This doesn't seem to work. So, I think that the generalized sum statistics approach of coin is not right for this case.
## We'd have to adapte edist() etc.. to the block randomized case. This alone might be a small article somewhere.
it2 <- independence_test(eY ~ newZ | bF,
  data = idat,
  distribution = approximate(nresample = 1000, parallel = "multicore", ncpus = 4),
  alternative = "greater"
)
size(it2, alpha = .05)

thepse <- matrix(NA, ncol = nsims, nrow = 2)
for (i in 1:nsims) {
  idat[, newZ := sample(Z), by = bF]
  idat[, eY := edisti(Ynull, newZ), by = bF]
  thepse[1, i] <- pvalue(independence_test(eY ~ newZ | bF,
    data = idat,
    distribution = approximate(nresample = 1000, parallel = "multicore", ncpus = 4)
  ))
  thepse[2, i] <- pvalue(independence_test(eY ~ newZ | bF,
    data = idat, alternative = "greater",
    distribution = approximate(nresample = 1000, parallel = "multicore", ncpus = 4)
  ))
}
apply(thepse, 1, function(x) {
  mean(x < .05)
})
