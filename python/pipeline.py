"""Loads the data and builds the objects of steps 2 and 3 (counterpart of R/pipeline.R).

Run from the repository root with python/ on the path:
    from pipeline import *
"""
import matplotlib.pyplot as plt  # noqa: F401 (re-exported)
import numpy as np  # noqa: F401
import pandas as pd
from dgxe import *  # noqa: F403

wx = pd.read_csv("data/narrabri_weather.csv")
traj = pd.read_csv("data/narrabri_biomass.csv.gz")
drv = make_drivers(wx)  # noqa: F405
Theta = make_library(drv)  # noqa: F405
tr = make_trajectories(traj, drv)  # noqa: F405
