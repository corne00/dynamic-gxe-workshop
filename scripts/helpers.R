## Helpers for the build scripts (not part of the code shown to the reader).

## n genotypes spread evenly over the ranking of mean final biomass.
spread_genotypes <- function(tr, n = 3) {
  last <- cbind(as.vector(row(tr$nend)), as.vector(col(tr$nend)),
                as.vector(tr$nend) + 1)
  final <- matrix(tr$Y[last], nrow(tr$nend), dimnames = dimnames(tr$nend))
  ranked <- names(sort(rowMeans(final)))
  ranked[round(seq(1, length(ranked), length.out = n))]
}
