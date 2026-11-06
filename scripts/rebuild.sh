#!/bin/sh
## Rebuild the site from scratch: precomputed results, figures, one render.
## _freeze is removed because Quarto does not notice changes in R/ and python/.
set -e
rm -rf _freeze precomputed/*.rds
Rscript scripts/precompute.R
Rscript scripts/concept_figures.R
quarto render
