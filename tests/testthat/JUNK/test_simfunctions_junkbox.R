##### OLD STUFF BELOW: Still should look at the errsimfn stuff and update it.



### Calculating error rates
res_biggrp <- norm_inc_det[, .(
  nodenum = paste(unique(nodenum_current), collapse = ","),
  blocks = paste(unique(bF), collapse = ","),
  nbiggrp = length(unique(biggrp)), # stri_count(biggrp,fixed=",")
  nblocks = length(unique(bF)),
  numtruenull = sum(ate_norm_inc == 0),
  anytruenull = as.numeric(any(ate_norm_inc == 0)),
  numtruetau = sum(ate_norm_inc != 0),
  anytruetau = as.numeric(any(ate_norm_inc != 0)),
  mntau = mean(ate_norm_inc),
  mintau = min(ate_norm_inc),
  maxtau = max(ate_norm_inc),
  meannb = mean(nb),
  minnb = min(nb),
  maxnb = max(nb),
  maxalpha = max(max_alpha),
  minalpha = min(max_alpha),
  hit = unique(hit),
  hit_grp = unique(hit_grp)
), by = biggrp][order(anytruetau, hit_grp), ]
res_biggrp
## All hits for true effects (no hits for groups where at least one of the blocks had a true effect)
with(res_biggrp, table(hit, truetau, exclude = c()))


res_hitgrp <- norm_inc_det[, .(
  nodenum = paste(unique(nodenum_current), collapse = ","),
  blocks = paste(unique(bF), collapse = ","),
  biggrps = abbreviate(paste(unique(biggrp), collapse = ",")),
  nbiggrp = length(unique(biggrp)), # stri_count(biggrp,fixed=",")
  nblocks = length(unique(bF)),
  hit = unique(hit),
  # p = min(pfinalb),
  # maxp=min(max_p),
  numtruenull = sum(ate_norm_inc == 0),
  anytruenull = as.numeric(any(ate_norm_inc == 0)),
  numtruetau = sum(ate_norm_inc != 0),
  anytruetau = as.numeric(any(ate_norm_inc != 0)),
  mntau = mean(ate_norm_inc),
  mintau = min(ate_norm_inc),
  maxtau = max(ate_norm_inc),
  meannb = mean(nb),
  minnb = min(nb),
  maxnb = max(nb),
  maxalpha = max(max_alpha),
  minalpha = min(max_alpha)
),
by = hit_grp
][order(anytruetau, hit_grp), ]
res_hitgrp

break


## Summarize errors across the splits / blocks
##
res3 <- res_biggrp[, .(
  detection_rate_fper = mean(as.numeric(any(p <= .05))), ## did we make any detections among all of the tests (each biggrp is a test)
  detections = as.numeric(any(hit)), ## did we make any detections among all of the tests (each biggrp is a test)
  fdp = sum(as.numeric(p <= .05) * (1 - truetau)) / max(1, sum(as.numeric(p <= .05))), ## false discovery proportion (false detections/rejections out of total rejections/detections)
  false_detect = mean(hit * (1 - truetau)), ## proportion hits when tau==0
  true_detect = mean(hit * truetau), ## proportion hits when tau!=0
  tdp = sum(as.numeric(p <= .05) * truetau) / max(1, sum(truetau)), ## true discovery proportion out of total true nonzero tau
  ngrps = .N
)]

## All hits for true effects (no hits for groups where at least one of the blocks had a true effect)
with(res_biggrp, table(hit, truetau, exclude = c()))

with(norm_inc_det, table(hit, ate_norm_inc == 0, exclude = c()))

res_biggrp[hit_grp == 8, ]
norm_inc_det[hit_grp == 8, ]



