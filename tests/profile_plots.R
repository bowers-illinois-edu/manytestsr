## Summarize and learn from the benchmarking exercise
### Overall, byunit_arma2_par  is best above about 20 and arma below
### If parallel is not available then overall the "arma" algorithmn wins below about N=200 and byunit_arma wins above

library(tidyverse)
library(bench)

load("dist_timings_small_osxlaptop.rda")
load("dist_timings_small_osxpro.rda")
load("dist_timings_small_linuxkeeling.rda")
load("dist_timings_large_osxlaptop.rda")
load("dist_timings_large_osxpro.rda")
load("dist_timings_large_linuxkeeling.rda")

## pdf(file="dist_timings_osxlaptop.pdf")
## ggplot2::autoplot(dist_timings_osxlaptop)
## dev.off()
#
osxprotime_small <- select(dist_timings_small_osxpro, c("N", "min", "median", "memory"))
osxprotime_small$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 19)
osxtime_small <- select(dist_timings_small_osxlaptop, c("N", "min", "median", "memory"))
osxtime_small$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 19)
linuxtime_small <- select(dist_timings_small_linuxkeeling, c("N", "min", "median", "memory"))
linuxtime_small$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 19)

osxprotime_large <- select(dist_timings_large_osxpro, c("N", "min", "median", "memory"))
osxprotime_large$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 11)
osxtime_large <- select(dist_timings_large_osxlaptop, c("N", "min", "median", "memory"))
osxtime_large$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 11)
linuxtime_large <- select(dist_timings_large_linuxkeeling, c("N", "min", "median", "memory"))
linuxtime_large$algo <- rep(c("main", "arma", "byunit_arma", "byunit_arma2", "byunit_arma2_par"), 11)

## If parallel is allowed, then anything above 20 should be parallelized using byunit_arma2_par, under 20 is arma
linuxtime_small %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxtime_small %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxprotime_small %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))

## If parallel is allowed, then anything above 1000 should be parallelized using byunit_arma2_par
linuxtime_large %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxtime_large %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxprotime_large %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))


## If parallel is not allowed:
### on linux arma until 100, main 100 to 1000 (or byunit_arma)
## on osx laptop arma until 400, main 400 to 1000 (or byunit_arma at 1000)
## on osx pro arma until about 200, main 200 to 500, byunit_arma 600+
linuxtime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxtime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))
osxprotime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  group_by(N) %>%
  mutate(newtime = (min + median) / 2) %>%
  filter(newtime == min(newtime))

linuxnonpar_small <- linuxtime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)
osxnonpar_small <- osxtime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)
osxprononpar_small <- osxprotime_small %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)

linuxnonpar_large <- linuxtime_large %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)
osxnonpar_large <- osxtime_large %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)
osxprononpar_large <- osxprotime_large %>%
  filter(algo != "byunit_arma2_par") %>%
  arrange(N, median) # %>%  print(n=100)

## devtools::install_github("thomasp85/patchwork")
library(patchwork)

glinuxnonpar_small1 <- ggplot(data = filter(linuxnonpar_small, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

glinuxnonpar_small2 <- ggplot(data = filter(linuxnonpar_small, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## Looks to me like byunit_arma until 1000 (we could do arma below 100 if we wanted)
glinuxnonpar_small1 + glinuxnonpar_small2


####

gosxnonpar_small1 <- ggplot(data = filter(osxnonpar_small, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxnonpar_small2 <- ggplot(data = filter(osxnonpar_small, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## Looks to me like arma until 600, and then main to 1000 (but could easily stop arma at 200)
gosxnonpar_small1 + gosxnonpar_small2

####

gosxprononpar_small1 <- ggplot(data = filter(osxprononpar_small, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxprononpar_small2 <- ggplot(data = filter(osxprononpar_small, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## Looks to me like arma until 200, main to 600, byunit_arma until 1000
gosxprononpar_small1 + gosxprononpar_small2
## Overall, looks like arma is best or at least very close to best below 200 and then byunit_arma is best (or very close to best) for the rest of the range.

##########################
### The larger vectors

glinuxnonpar_large1 <- ggplot(data = filter(linuxnonpar_large, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

glinuxnonpar_large2 <- ggplot(data = filter(linuxnonpar_large, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## byunit_arma is best in this range with byunit_arma2 a close second
glinuxnonpar_large1 + glinuxnonpar_large2

####

gosxnonpar_large1 <- ggplot(data = filter(osxnonpar_large, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxnonpar_large2 <- ggplot(data = filter(osxnonpar_large, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## byunit_arma is best in this range with byunit_arma2 a very very close second
gosxnonpar_large1 + gosxnonpar_large2

####

gosxprononpar_large1 <- ggplot(data = filter(osxprononpar_large, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxprononpar_large2 <- ggplot(data = filter(osxprononpar_large, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## byunit_arma is best in this range with byunit_arma2 a close second
gosxprononpar_large1 + gosxprononpar_large2


##########################
###  All algos

dist_timings_small_linuxkeeling$algo <- names(dist_timings_small_linuxkeeling$expression)
dist_timings_small_osxlaptop$algo <- names(dist_timings_small_osxlaptop$expression)
dist_timings_small_osxpro$algo <- names(dist_timings_small_osxpro$expression)
dist_timings_large_linuxkeeling$algo <- names(dist_timings_large_linuxkeeling$expression)
dist_timings_large_osxlaptop$algo <- names(dist_timings_large_osxlaptop$expression)
dist_timings_large_osxpro$algo <- names(dist_timings_large_osxpro$expression)

###
glinux_small1 <- ggplot(data = filter(dist_timings_small_linuxkeeling, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

glinux_small2 <- ggplot(data = filter(dist_timings_small_linuxkeeling, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## arma wins but only under about 20
glinux_small1 + glinux_small2

####

gosxlaptop_small1 <- ggplot(data = filter(dist_timings_small_osxlaptop, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxlaptop_small2 <- ggplot(data = filter(dist_timings_small_osxlaptop, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## same arma but only under 20
gosxlaptop_small1 + gosxlaptop_small2


####

gosxpro_small1 <- ggplot(data = filter(dist_timings_small_osxpro, N <= 200), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxpro_small2 <- ggplot(data = filter(dist_timings_small_osxpro, N >= 200), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

## same
gosxpro_small1 + gosxpro_small2


###
glinux_large1 <- ggplot(data = filter(dist_timings_large_linuxkeeling, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

glinux_large2 <- ggplot(data = filter(dist_timings_large_linuxkeeling, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

glinux_large1 + glinux_large2

####

gosxlaptop_large1 <- ggplot(data = filter(dist_timings_large_osxlaptop, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxlaptop_large2 <- ggplot(data = filter(dist_timings_large_osxlaptop, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

gosxlaptop_large1 + gosxlaptop_large2


####

gosxpro_large1 <- ggplot(data = filter(dist_timings_large_osxpro, N <= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line(show.legend = FALSE) +
  theme_bw()

gosxpro_large2 <- ggplot(data = filter(dist_timings_large_osxpro, N >= 5000), aes(N, median, group = algo, color = algo)) +
  geom_line() +
  theme(legend.position = c(400, .01), legend.background = element_blank()) +
  theme_bw()

gosxpro_large1 + gosxpro_large2
