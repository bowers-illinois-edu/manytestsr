

library(bench)
library(tidyverse)
load("bench_press_results_dist_timings.rda")

summary(dist_timings,filter_gc=FALSE)

ggplot2::autoplot(filter(dist_timings,n_gc==0))

