## All precomputed results in one pass: the cross-validation and ALE results
## that the pages read from precomputed/, and the data for the interactive
## figures (data/interactive/*.json). Run from the repository root, before
## rendering: Rscript scripts/precompute.R   (or scripts/rebuild.sh for all).
source("R/pipeline.R")
source("R/fast.R")
suppressMessages(library(iml))
dir.create("precomputed", showWarnings = FALSE)
store <- function(file, x) saveRDS(x, file.path("precomputed", file))
cores <- parallel::detectCores()

## Held-out predictions of the sparse model, and the method comparison (step 5).
## OLS, ridge and STLSQ use the moment-based fits (same values as fit_ols,
## fit_ridge and fit_stlsq); the STLSQ predictions serve both results.
moment_fit <- list(OLS = list(mpath_ols, 0), ridge = list(mpath_ridge, 1e-3),
                   STLSQ = list(mpath_stlsq, 0.05))
P_moment <- lapply(moment_fit, function(m) {
  cv_moments(Theta, tr, m[[1]], m[[2]], what = "P")[[1]]
})
store("cv_stlsq.rds", P_moment$STLSQ)
store("cv_methods.rds", t(sapply(names(fit_methods), function(name) {
  P <- if (name %in% names(P_moment)) P_moment[[name]]
       else cross_validate_par(Theta, tr, fit_methods[[name]])
  evaluate(P, tr)
})))

## Error against the number of terms, with products of pairs (step 5)
Theta2 <- make_library(drv, pairs = TRUE)
ks <- c(4, 8, 12, 16, 24, 36, 60, 108)
store("rmse_vs_k.rds", data.frame(k = ks, rmse = cv_rmse_path(
  Theta2, tr, path_subset("backward"), ks, tr$genotypes[seq(1, 199, by = 20)])))

## ALE curves of every genotype and driver (step 6)
coefs <- sapply(tr$genotypes, function(g) {
  d <- make_design(Theta, tr, g)
  fit_stlsq(d$X, d$dy)
})
Z <- as.data.frame(observed_points(drv, tr))
predictor <- function(coef, pairs = FALSE) {
  rate <- function(model, newdata) {
    L <- library_matrix(as.matrix(newdata), pairs = pairs)
    as.vector(cbind(1, L) %*% model)
  }
  Predictor$new(coef, data = Z, predict.function = rate)
}
natural <- function(z, j) z * drv$scale[j] + drv$center[j]
ale_genotype <- function(g, j) {
  effect <- FeatureEffect$new(predictor(coefs[, g]), drv$drivers[j],
                              method = "ale", grid.size = 24)$results
  data.frame(genotype = g, driver = factor(drv$drivers[j], drv$drivers),
             x = natural(effect[[drv$drivers[j]]], j), effect = effect$.value)
}
cells <- expand.grid(g = tr$genotypes, j = seq_along(drv$drivers),
                     stringsAsFactors = FALSE)
store("ale_genotypes.rds", do.call(rbind, parallel::mclapply(
  seq_len(nrow(cells)), function(i) ale_genotype(cells$g[i], cells$j[i]),
  mc.cores = cores)))

source("scripts/precompute_interactive.R")
