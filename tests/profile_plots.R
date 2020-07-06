## Summarize and learn from the benchmarking exercise

library(tidyverse)
library(bench)

load("dist_timings_small_osxlaptop.rda")
load("dist_timings_small_osxpro.rda")
load("dist_timings_small_linuxkeeling.rda")
load("dist_timings_large_osxlaptop.rda")
load("dist_timings_large_osxpro.rda")
load("dist_timings_large_linuxkeeling.rda")

print(dist_timings_osxlaptop, n=100)
## Something about fast_dists_and_trans_by_unit is just destroying memory. The byunit_arma is the right approach for big stuff.

pdf(file="dist_timings_osxlaptop.pdf")
ggplot2::autoplot(dist_timings_osxlaptop)
dev.off()


pdf(file="dist_timings_linuxkeeling.pdf")
ggplot2::autoplot(dist_timings_linuxkeeling)
dev.off()

print(dist_timings_linuxkeeling, n=100)
print(dist_timings_osxlaptop, n=100)


osxtime <- select(dist_timings_osxlaptop,c("N","min","median","memory"))
osxtime$algo <- rep(c("main","arma","byunit_arma","byunit_arma2","byunit_arma2_par"),8)

linuxtime <- select(dist_timings_linuxkeeling,c("N","min","median","memory"))
linuxtime$algo <- rep(c("main","arma","byunit_arma","byunit_arma2","byunit_arma2_par"),8)


## If parallel is allowed, then anything above 100 should be parallelized using byunit_arma2_par
## Under 100 use byunit_arma

linuxtime %>% group_by(N) %>% mutate(newtime = (min+median)/2) %>% filter(newtime==min(newtime))
osxtime %>% group_by(N) %>% mutate(newtime = (min+median)/2) %>% filter(newtime==min(newtime))


## With out parallelization, then it is mostly byunit_arma. (sometimes tied with arma). Even in when N=500, main is basically tied with byunit_arma.

linuxtime %>% filter(algo != "byunit_arma2_par") %>% group_by(N) %>% mutate(newtime = (min+median)/2) %>% filter(newtime==min(newtime))


print(linuxtime %>% filter(algo != "byunit_arma2_par"),n=100)

## Seems like the osx laptop prefers byunit_arma2 but only very slightly.

osxtime %>% filter(algo != "byunit_arma2_par") %>% group_by(N) %>% mutate(newtime = (min+median)/2) %>% filter(newtime==min(newtime))

print(osxtime %>% filter(algo != "byunit_arma2_par"),n=100)

## So: use arma N<50, and when parallel flag is set in pIndepDist byunit_arma2_par for N>50 or byunit_arma.
## parallel seems to be winning starting around N=20.

## But when we do a bit more benchmarking


