# Leave-one-year-out validation for a few genotypes and all fitting methods.
source("R/functions.R")
weather <- read.csv("data/narrabri_weather.csv")
traj <- read.csv("data/narrabri_biomass.csv.gz")
drv <- make_drivers(weather)
Theta <- make_library(drv)
tr <- make_trajectories(traj, drv)
genotypes <- tr$genotypes[seq(1, 199, by = 10)]
for (name in names(fit_methods)) {
  P <- cross_validate(Theta, tr, fit_methods[[name]], genotypes)
  cat(sprintf("%-22s", name), round(evaluate(P, tr), 3), "\n")
}
