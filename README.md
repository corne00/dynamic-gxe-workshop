# Dynamic G×E workshop

Hands-on workshop (AAGI 2026 Symposium): learning growth equations from crop-model data with sparse regression, using base R only.

Site: step-by-step pages (method, data, design matrix, fitting, validation, interpretation, genotypes, limitations, extensions) with copyable code and pre-run output.

## Run it yourself

```r
source("R/functions.R")
wx   <- read.csv("data/narrabri_weather.csv")
traj <- read.csv("data/narrabri_biomass.csv.gz")
```

`scripts/parity_check.R` runs the whole pipeline in a few seconds. Base R only; the pages additionally use `ggplot2` (and `splines` for one example).

## Build the site

Requires [Quarto](https://quarto.org) and R with `ggplot2`.

```sh
quarto render      # writes docs/ ; computed outputs are cached in _freeze/ and precomputed/
```

Publish with GitHub Pages: Settings → Pages → Deploy from a branch → `main`, folder `/docs`.

## Layout

- `*.qmd`: the pages · `R/functions.R`: all workshop functions · `R/pipeline.R`: steps 2–4 in one script
- `data/`: Narrabri weather and APSIM biomass/LAI trajectories · `scripts/`: data export and parity check
- `images/`: static figures from the research code (page 9 only)

## Data

The data are APSIM simulation output (199 genotypes × 31 seasons, Narrabri). The source and licence of the underlying simulations still need to be confirmed before reuse; no licence is granted for the data here.
