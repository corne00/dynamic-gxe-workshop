"""Leave-one-season-out validation for a few genotypes and all fitting
methods (counterpart of parity_check.R). Run from the repository root with
PYTHONPATH=python."""
import pandas as pd
from dgxe import *

wx = pd.read_csv("data/narrabri_weather.csv")
traj = pd.read_csv("data/narrabri_biomass.csv.gz")
drv = make_drivers(wx)
Theta = make_library(drv)
tr = make_trajectories(traj, drv)

genotypes = tr.genotypes[::10]
for name, fit in fit_methods.items():
    P = cross_validate(Theta, tr, fit, genotypes)
    print(f"{name:22s}", *(f"{v:.3f}" for v in evaluate(P, tr).values()))
