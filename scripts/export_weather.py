"""Narrabri weather table for the R workshop.
Given (simulated weather inputs): Tmin, Tmax, Rain, RAD, VP, DL.
Constructed from them: TT (accumulated thermal time) and RAW (relative available water).
Usage: python export_weather.py weather_daily_by_environment.csv model_drivers_daily.csv narrabri_weather.csv"""
import sys
import pandas as pd

weather, drivers, out = sys.argv[1:4]
w = pd.read_csv(weather).query("loc == 'Narrabri'")
d = pd.read_csv(drivers).query("loc == 'Narrabri'")
d = d[["Env", "das", "fao56_relative_available_water"]]
x = w.merge(d.assign(das=d.das.astype(float)), on=["Env", "das"])
x = x.rename(columns={"MinT": "Tmin", "MaxT": "Tmax", "Radn": "RAD", "vp": "VP", "day_length": "DL",
                      "temperature_integral": "TT", "fao56_relative_available_water": "RAW"})
x["das"] = x.das.astype(int)
x = x.sort_values(["year", "das"])[["year", "das", "Tmin", "Tmax", "Rain", "RAD", "VP", "DL", "TT", "RAW"]]
x.round({"TT": 2, "RAW": 4}).to_csv(out, index=False)
print(x.shape)
