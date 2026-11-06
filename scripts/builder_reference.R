## Reference results of the R functions for a few configurations, used to check
## js/builder.js (run: Rscript scripts/builder_reference.R; node scripts/builder_check.js).
library(jsonlite)
source("R/functions.R")
wx <- read.csv("data/narrabri_weather.csv"); traj <- read.csv("data/narrabri_biomass.csv.gz")
all8 <- c("Tmin", "Tmax", "Rain", "RAD", "VP", "DL", "TT", "RAW")
configs <- list(
  list(drivers = all8, degree = 3, pairs = FALSE, method = "OLS", value = NULL),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "ridge", value = 1e-3),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "lasso", value = 3e-4),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "elastic net", value = 6e-4),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "STLSQ", value = 0.05),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "backward", value = 8),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "best subset", value = 5),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "STLSQ", value = 0.03),
  list(drivers = c("TT", "DL", "RAW"), degree = 2, pairs = TRUE, method = "backward", value = 10),
  list(drivers = c("Tmax", "RAD", "RAW"), degree = 1, pairs = FALSE, method = "lasso", value = 1e-3))
fitter <- function(m, v) switch(m,
  OLS = fit_ols, ridge = function(X, dy) fit_ridge(X, dy, v),
  lasso = function(X, dy) fit_lasso(X, dy, v), `elastic net` = function(X, dy) fit_enet(X, dy, v),
  STLSQ = function(X, dy) fit_stlsq(X, dy, v),
  backward = function(X, dy) fit_subset(X, dy, v, "backward"),
  `best subset` = function(X, dy) fit_subset(X, dy, v))
genotypes <- jsonlite::read_json("data/interactive/builder.json")$genotypes |> unlist()
out <- lapply(configs, function(cf) {
  drv <- make_drivers(wx, cf$drivers); Theta <- make_library(drv, cf$degree, cf$pairs)
  tr <- make_trajectories(traj, drv)
  g <- genotypes[c(1, 5)]
  res <- lapply(g, function(gg) {
    f <- fitter(cf$method, cf$value)
    P <- cross_validate(Theta, tr, f, gg)
    d <- make_design(Theta, tr, gg); b <- f(d$X, d$dy)
    list(genotype = gg, rmse = unname(evaluate(P, tr)["trajectory_rmse"]),
         maturity = unname(evaluate(P, tr)["maturity_rmse"]), nterms = sum(b[-1] != 0))
  })
  c(cf, list(results = res))
})
write_json(out, "scripts/builder_reference.json", auto_unbox = TRUE, digits = NA, null = "null")
