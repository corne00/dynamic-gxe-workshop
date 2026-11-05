# Dynamic G×E workshop

Hands-on workshop (AAGI 2026 Symposium): learning growth equations from crop-model data with sparse regression, in R and Python.

Site: step-by-step pages (method, data, design matrix, fitting, validation, interpretation, genotypes, limitations, extensions) with copyable code and pre-run output.

## Run it yourself

```r
source("R/functions.R")
wx   <- read.csv("data/narrabri_weather.csv")
traj <- read.csv("data/narrabri_biomass.csv.gz")
```

Python: `PYTHONPATH=python`, then `from dgxe import *`. `scripts/parity_check.R` and `scripts/parity_check.py` cross-validate all fitting methods on a few genotypes and give the same errors. R packages: `glmnet`, `leaps`, `iml`, `ggplot2` (and `splines` for one example). Python packages: `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `PyALE`.

## Build the site

Requires [Quarto](https://quarto.org) and R with `ggplot2`, `glmnet`, `leaps` and `iml`.

```sh
quarto render      # writes docs/ ; computed outputs are cached in _freeze/ and precomputed/
```

Publish with GitHub Pages: Settings → Pages → Deploy from a branch → `main`, folder `/docs`.

## Interactive figures

Pages 1 and 4–7 contain sliders and menus (Observable JS, with `d3` and Observable Plot in `js/`, no external requests). They read results computed beforehand by `scripts/precompute_interactive.R` (a few minutes) and `scripts/precompute_weather.R`; both write `data/interactive/*.json`. Page 1 has a toy of the chain environment, response and biomass, driven by the observed temperature of a Narrabri season.

## Layout

- `*.qmd`: the pages · `R/functions.R` and `python/dgxe.py`: all workshop functions · `R/pipeline.R`: steps 2–3 in one script
- `data/`: Narrabri weather and APSIM biomass/LAI trajectories · `scripts/`: data export and parity check
- `images/`: static figures from the research code (page 9 only)

## Data

The data are APSIM simulation output (199 genotypes × 31 seasons, Narrabri). The source and licence of the underlying simulations still need to be confirmed before reuse; no licence is granted for the data here.
