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
  list(drivers = c("TT", "DL", "RAW"), degree = 2, pairs = TRUE, method = "backward", value = 7),
  list(drivers = c("Tmax", "RAD", "RAW"), degree = 1, pairs = FALSE, method = "lasso", value = 1e-3),
  ## measured trait on the right-hand side, and LAI as the modelled state
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "STLSQ", value = 0.03,
       target = "biomass", trait = "lai"),
  list(drivers = all8, degree = 3, pairs = FALSE, method = "lasso", value = 3e-4,
       target = "biomass", trait = "lai"),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "ridge", value = 1e-3,
       target = "lai"),
  list(drivers = c("TT", "RAD", "DL"), degree = 2, pairs = TRUE, method = "backward", value = 8,
       target = "lai", trait = "biomass"),
  list(drivers = c("TT", "DL", "RAW"), degree = 1, pairs = FALSE, method = "OLS", value = NULL,
       target = "biomass", trait = "lai", extras = c(canopy = "lai * RAD", wet = "roll(lai, 5) - RAW")))
## coupled models and state feedback: states, feedback
configs <- c(configs, list(
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "STLSQ", value = 0.03,
       states = c("biomass", "lai"), feedback = FALSE),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "STLSQ", value = 0.03,
       states = c("biomass", "lai"), feedback = TRUE),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "STLSQ", value = 0.03,
       states = "biomass", feedback = TRUE),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "ridge", value = 1e-3,
       states = "lai", feedback = TRUE),
  list(drivers = c("TT", "RAD"), degree = 2, pairs = TRUE, method = "backward", value = 8,
       states = c("biomass", "lai"), feedback = TRUE),
  list(drivers = c("TT", "DL", "RAW"), degree = 2, pairs = FALSE, method = "lasso", value = 1e-3,
       states = "biomass", feedback = TRUE, trait = "lai"),
  list(drivers = c("TT", "DL", "RAW", "RAD"), degree = 2, pairs = FALSE, method = "lasso", value = 3e-3,
       states = c("biomass", "lai"), feedback = TRUE),
  list(drivers = c("TT", "RAD", "RAW"), degree = 2, pairs = TRUE, method = "ridge", value = 0.1,
       states = c("biomass", "lai"), feedback = TRUE),
  list(drivers = c("TT", "DL", "RAW"), degree = 3, pairs = TRUE, method = "backward", value = 10),
  list(drivers = c("TT", "RAD", "RAW"), degree = 3, pairs = TRUE, method = "ridge", value = 0.1,
       states = c("biomass", "lai"), feedback = TRUE)))
fitter <- function(m, v) switch(m,
  OLS = fit_ols, ridge = function(X, dy) fit_ridge(X, dy, v),
  lasso = function(X, dy) fit_lasso(X, dy, v), `elastic net` = function(X, dy) fit_enet(X, dy, v),
  STLSQ = function(X, dy) fit_stlsq(X, dy, v),
  backward = function(X, dy) fit_subset(X, dy, v, "backward"),
  `best subset` = function(X, dy) fit_subset(X, dy, v))
genotypes <- jsonlite::read_json("data/interactive/builder.json")$genotypes |> unlist()
out <- lapply(configs, function(cf) {
  target <- if (is.null(cf$target)) "biomass" else cf$target
  g <- genotypes[c(1, 5)]
  if (!is.null(cf$states)) {
    res <- lapply(g, function(gg) {
      wg <- if (is.null(cf$trait)) wx else with_traits(wx, traj, gg, cf$trait)
      all <- c(cf$drivers, cf$trait)
      drv <- make_drivers(wg, all); Theta <- make_library(drv, cf$degree, cf$pairs)
      trs <- setNames(lapply(cf$states, function(st) make_trajectories(traj, drv, st)), cf$states)
      f <- fitter(cf$method, cf$value)
      P <- cross_validate_states(Theta, trs, f, gg, all, cf$degree, cf$pairs, cf$feedback)
      ev <- lapply(seq_along(trs), function(k) evaluate(P[[k]], trs[[k]]))
      list(genotype = gg, rmse = sapply(ev, `[[`, "trajectory_rmse"),
           maturity = sapply(ev, `[[`, "maturity_rmse"), diverged = attr(P, "diverged"))
    })
    return(c(cf, list(results = res)))
  }
  res <- lapply(g, function(gg) {
    wg <- if (is.null(cf$trait)) wx else with_traits(wx, traj, gg, cf$trait)
    drv <- make_drivers(wg, c(cf$drivers, cf$trait)); Theta <- make_library(drv, cf$degree, cf$pairs)
    if (!is.null(cf$extras)) Theta <- bind_terms(Theta, extra_terms(wg, cf$extras))
    tr <- make_trajectories(traj, drv, target)
    f <- fitter(cf$method, cf$value)
    P <- cross_validate(Theta, tr, f, gg)
    d <- make_design(Theta, tr, gg); b <- f(d$X, d$dy)
    list(genotype = gg, rmse = unname(evaluate(P, tr)["trajectory_rmse"]),
         maturity = unname(evaluate(P, tr)["maturity_rmse"]), nterms = sum(b[-1] != 0))
  })
  c(cf, list(results = res))
})
write_json(out, "scripts/builder_reference.json", auto_unbox = TRUE, digits = NA, null = "null")
