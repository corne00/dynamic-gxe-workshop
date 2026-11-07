#!/bin/sh
## Rebuild the site. Each stage runs only when its inputs changed since the last
## run (checksums in .build/); "scripts/rebuild.sh --force" recomputes everything.
##   1. precomputed results (precomputed/, data/interactive/)  [slowest]
##   2. concept figures (images/)
##   3. render; _freeze is removed when R/ or python/ changed, because Quarto
##      does not notice those, and otherwise reuses the unchanged pages.
set -e
cd "$(dirname "$0")/.."
[ "$1" = "--force" ] && rm -rf .build _freeze precomputed/*.rds
mkdir -p .build
sum_of() { cat "$@" | md5sum | cut -d' ' -f1; }
stale() { [ "$(cat ".build/$1" 2>/dev/null)" != "$2" ]; }
data="data/narrabri_weather.csv data/narrabri_biomass.csv.gz"
code="R/functions.R R/pipeline.R R/fast.R"

s1=$(sum_of $code $data scripts/precompute.R scripts/precompute_interactive.R)
if stale precompute "$s1" || [ ! -f precomputed/cv_methods.rds ]; then
  Rscript scripts/precompute.R
  echo "$s1" > .build/precompute
fi

s2=$(sum_of $code $data scripts/concept_figures.R R/theme_dgxe.R)
if stale figures "$s2"; then
  Rscript scripts/concept_figures.R
  echo "$s2" > .build/figures
fi

s3=$(sum_of R/*.R python/*.py $data)
if stale render "$s3"; then
  rm -rf _freeze
fi
quarto render
echo "$s3" > .build/render
