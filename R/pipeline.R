## Steps 2-4 in one go (used by the later pages so that each page runs on its own).
source("R/functions.R")
wx    <- read.csv("data/narrabri_weather.csv")
traj  <- read.csv("data/narrabri_biomass.csv.gz")
drv   <- make_drivers(wx)
Theta <- make_library(drv)
tr    <- make_trajectories(traj, drv)
gram  <- gram_stats(Theta, tr)
