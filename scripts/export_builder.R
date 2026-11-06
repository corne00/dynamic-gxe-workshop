## Data for the "build your own model" page (data/interactive/builder.json):
## standardised drivers of all years and the biomass of a few genotypes.
library(jsonlite)
source("R/pipeline.R")

genotypes <- spread_genotypes(tr, 8)
pick <- match(genotypes, tr$genotypes)
write_json(list(
  drivers = drv$drivers, center = drv$center, scale = drv$scale,
  years = tr$years, das0 = tr$das0, genotypes = genotypes,
  U = lapply(seq_along(tr$years), function(y) {
    lapply(seq_len(dim(drv$U)[2]), function(d) round(drv$U[y, d, ], 4))
  }),
  nend = lapply(pick, function(g) tr$nend[g, ]),
  Y = lapply(pick, function(g) lapply(seq_along(tr$years), function(y) {
    round(tr$Y[g, y, 1:(tr$nend[g, y] + 1)], 4)
  }))),
  "data/interactive/builder.json", auto_unbox = TRUE, digits = NA)
