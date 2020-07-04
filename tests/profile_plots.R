

library(tidyverse)
library(bench)

load("dist_timings_osxlaptop.rda")
load("dist_timings_linuxkeeling.rda")

print(dist_timings_osxlaptop, n=100)
## Something about fast_dists_and_trans_by_unit is just destroying memory. The byunit_arma is the right approach for big stuff.

pdf(file="dist_timings_osxlaptop.pdf")
ggplot2::autoplot(dist_timings_osxlaptop)
dev.off()

print(dist_timings_linuxkeeling, n=100)


## Looks like more or less the following from both machines:
## Under about 500 arma is best
## Between about 1000 and about 50000 main is best
## Above about 5000 byunit_arma is best









