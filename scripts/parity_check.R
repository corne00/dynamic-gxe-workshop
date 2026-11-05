# Runs the whole pipeline once (leave-one-year-out, shared support of k terms + genotype deviations)
source("R/functions.R")
t0 <- proc.time()[3]
wx   <- read.csv("data/narrabri_weather.csv")
traj <- read.csv(gzfile("data/narrabri_biomass.csv.gz"))
drv   <- make_drivers(wx); Theta <- make_library(drv); tr <- make_trajectories(traj, drv)
gram  <- gram_stats(Theta, tr); cat(sprintf("setup %.1fs; Theta %s; Y %s\n", proc.time()[3] - t0, paste(dim(Theta), collapse="x"), paste(dim(tr$Y), collapse="x")))
for (k in c(12)) {
  t1 <- proc.time()[3]; cv <- cv_blocked(Theta, tr, gram, fit_ssr, k = k)
  m <- evaluate(cv$pred, tr, fixed_day = 99)    # array index 99 = day index 98 = day 113 after sowing
  cat(sprintf("k=%d  trajectory RMSE %.3f   day-113 RMSE %.3f   maturity RMSE %.3f   (%.1fs)\n", k, m[1], m[2], m[3], proc.time()[3] - t1))
}
