## Loads the data and builds the objects of steps 2 and 3.
source("R/functions.R")
weather <- read.csv("data/narrabri_weather.csv")
traj <- read.csv("data/narrabri_biomass.csv.gz")
drv <- make_drivers(weather)
Theta <- make_library(drv)
tr <- make_trajectories(traj, drv)
