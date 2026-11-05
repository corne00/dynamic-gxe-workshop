## =============================================================================
##  Sparse identification of crop growth: functions used in the workshop
##
##    weather --make_drivers--> drv --make_library--> Theta   (candidate terms)
##    biomass --make_trajectories--> tr
##    Theta, tr --make_design--> X, dy                (regression of increments)
##    fit_*(X, dy) --> coefficients --simulate--> prediction --evaluate--> RMSE
##
##  Base R, except for glmnet and leaps in some of the fitting methods.
## =============================================================================

## ---- Drivers: array U[season, day, driver], standardised -------------------
make_drivers <- function(wx, drivers = c("TT", "DL", "RAW", "RAD")) {
  years <- sort(unique(wx$year))
  das0 <- min(wx$das)
  nd <- max(wx$das) - das0 + 1
  U <- array(NA_real_, c(length(years), nd, length(drivers)),
             dimnames = list(years, NULL, drivers))
  for (j in seq_along(drivers)) {
    i <- cbind(match(wx$year, years), wx$das - das0 + 1, j)
    U[i] <- wx[[drivers[j]]]
  }
  # seasons end on different days: pad with the last value (never fitted)
  for (y in seq_along(years)) {
    for (j in seq_along(drivers)) {
      u <- U[y, , j]
      U[y, is.na(u), j] <- tail(u[!is.na(u)], 1)
    }
  }
  center <- apply(U, 3, mean)
  scale <- apply(U, 3, function(x) sqrt(mean((x - mean(x))^2)))
  list(U = sweep(sweep(U, 3, center), 3, scale, "/"), years = years,
       das0 = das0, drivers = drivers, center = center, scale = scale)
}

## ---- Library of candidate terms ---------------------------------------------
## Z: n x m matrix of standardised drivers -> n x K matrix of terms. Powers
## 1..degree of every driver; pairs = TRUE adds products of pairs of drivers.
## The constant term is the intercept of the regression and is not included.
library_matrix <- function(Z, degree = 3, pairs = FALSE) {
  term <- function(x, p) if (p == 1) x else paste0(x, "^", p)
  cols <- list()
  for (d in colnames(Z)) {
    for (p in 1:degree) cols[[term(d, p)]] <- Z[, d]^p
  }
  if (pairs) {
    for (ij in combn(colnames(Z), 2, simplify = FALSE)) {
      for (p in 1:2) {
        for (q in 1:2) {
          name <- paste0(term(ij[1], p), "*", term(ij[2], q))
          cols[[name]] <- Z[, ij[1]]^p * Z[, ij[2]]^q
        }
      }
    }
  }
  do.call(cbind, cols)
}

## Theta[season, day, term]
make_library <- function(drv, degree = 3, pairs = FALSE) {
  d <- dim(drv$U)
  Z <- matrix(drv$U, ncol = d[3], dimnames = list(NULL, drv$drivers))
  L <- library_matrix(Z, degree, pairs)
  array(L, c(d[1:2], ncol(L)), dimnames = list(NULL, NULL, colnames(L)))
}

## ---- Trajectories: array Y[genotype, season, day], NA after maturity --------
## nend[g, y] is the number of daily steps of trajectory (g, y).
make_trajectories <- function(traj, drv, var = "biomass") {
  gen <- sort(unique(traj$genotype))
  Y <- array(NA_real_, c(length(gen), length(drv$years), dim(drv$U)[2]),
             dimnames = list(gen, drv$years, NULL))
  i <- cbind(match(traj$genotype, gen), match(traj$year, drv$years),
             traj$das - drv$das0 + 1)
  Y[i] <- traj[[var]]
  nend <- apply(Y, c(1, 2), function(v) sum(!is.na(v)) - 1)
  list(Y = Y, nend = nend, genotypes = gen, years = drv$years,
       das0 = drv$das0)
}

## ---- Regression problem for one genotype ------------------------------------
## Rows: all days of all seasons. X: library at day t. dy: y(t + 1) - y(t).
make_design <- function(Theta, tr, genotype) {
  g <- match(genotype, tr$genotypes)
  n <- tr$nend[g, ]
  seasons <- seq_along(tr$years)
  X <- lapply(seasons, function(y) matrix(Theta[y, 1:n[y], ], n[y]))
  dy <- lapply(seasons, function(y) diff(tr$Y[g, y, 1:(n[y] + 1)]))
  X <- do.call(rbind, X)
  colnames(X) <- dimnames(Theta)[[3]]
  list(X = X, dy = unlist(dy), season = rep(seasons, n))
}

## ---- Fitting methods: fit(X, dy) -> coefficients (intercept first) ----------
named <- function(b, X) setNames(unname(b), c("1", colnames(X)))

fit_ols <- function(X, dy) {
  lm.fit(cbind(`1` = 1, X), dy)$coefficients
}

## Ridge (alpha = 0), lasso (alpha = 1) and elastic net by glmnet.
fit_glmnet <- function(X, dy, alpha, lambda) {
  fit <- glmnet::glmnet(X, dy, alpha = alpha, lambda = lambda)
  named(as.vector(coef(fit)), X)
}
fit_ridge <- function(X, dy, lambda = 1e-3) fit_glmnet(X, dy, 0, lambda)
fit_lasso <- function(X, dy, lambda = 3e-4) fit_glmnet(X, dy, 1, lambda)
fit_enet <- function(X, dy, lambda = 6e-4) fit_glmnet(X, dy, 0.5, lambda)

## Best subset (exhaustive) or backward elimination with k terms, by leaps.
fit_subset <- function(X, dy, k = 8, method = "exhaustive") {
  b <- coef(leaps::regsubsets(X, dy, nvmax = k, method = method), k)
  full <- setNames(numeric(ncol(X) + 1), c("1", colnames(X)))
  full[c("1", names(b)[-1])] <- b
  full
}

## Sequential thresholded least squares (the original SINDy algorithm): refit
## on the active terms and remove those whose scaled coefficient is below the
## threshold, until no term is removed.
fit_stlsq <- function(X, dy, threshold = 0.05) {
  size <- apply(X, 2, sd) / sd(dy)         # puts all coefficients on one scale
  active <- colnames(X)
  repeat {
    b <- lm.fit(cbind(`1` = 1, X[, active, drop = FALSE]), dy)$coefficients
    keep <- abs(b[active]) * size[active] >= threshold
    if (all(keep) || !any(keep)) break
    active <- active[keep]
  }
  full <- setNames(numeric(ncol(X) + 1), c("1", colnames(X)))
  full[names(b)] <- b
  full
}

## Methods compared in the workshop (hyperparameters fixed for simplicity).
fit_methods <- list(
  OLS = fit_ols,
  ridge = fit_ridge,
  lasso = fit_lasso,
  `elastic net` = fit_enet,
  `best subset (k = 8)` = function(X, dy) fit_subset(X, dy, k = 8),
  `backward (k = 8)` = function(X, dy) fit_subset(X, dy, 8, "backward"),
  STLSQ = fit_stlsq)

## ---- Simulation: for a weather-only library the solution is a cumulative sum
## coef: intercept and term coefficients; Theta_y: day x term (one season).
simulate <- function(coef, Theta_y, y0) {
  rate <- cbind(1, Theta_y[-nrow(Theta_y), ]) %*% coef
  y0 + c(0, cumsum(rate))
}

## Genotypes with the lowest, median and highest mean final biomass.
spread_genotypes <- function(tr) {
  last <- cbind(as.vector(row(tr$nend)), as.vector(col(tr$nend)),
                as.vector(tr$nend) + 1)
  final <- matrix(tr$Y[last], nrow(tr$nend), dimnames = dimnames(tr$nend))
  ranked <- names(sort(rowMeans(final)))
  ranked[c(1, ceiling(length(ranked) / 2), length(ranked))]
}

## ---- Validation: leave-one-season-out ---------------------------------------
## Each season is predicted by a model fitted to the other seasons of the same
## genotype, starting from the mean initial biomass of those seasons.
## Returns P[genotype, season, day].
cross_validate <- function(Theta, tr, fit, genotypes = tr$genotypes) {
  P <- array(NA_real_, dim(tr$Y))
  for (g in match(genotypes, tr$genotypes)) {
    d <- make_design(Theta, tr, tr$genotypes[g])
    for (v in seq_along(tr$years)) {
      train <- d$season != v
      coef <- fit(d$X[train, ], d$dy[train])
      P[g, v, ] <- simulate(coef, Theta[v, , ], mean(tr$Y[g, -v, 1]))
    }
  }
  P
}

## RMSE over whole trajectories, on a fixed day (array index; 99 is day 113
## after sowing) and at maturity. Genotypes that were not predicted are ignored.
evaluate <- function(P, tr, fixed_day = 99) {
  err <- P - tr$Y
  err[, , 1] <- NA                                # the initial value is given
  last <- cbind(as.vector(row(tr$nend)), as.vector(col(tr$nend)),
                as.vector(tr$nend) + 1)
  rmse <- function(e) sqrt(mean(e^2, na.rm = TRUE))
  c(trajectory_rmse = rmse(err), fixed_day_rmse = rmse(err[, , fixed_day]),
    maturity_rmse = rmse(err[last]))
}

## RMSE on each day, over all genotype-season pairs still growing on that day.
rmse_by_day <- function(P, tr) {
  sqrt(apply((P - tr$Y)^2, 3, mean, na.rm = TRUE))[-1]
}

## Baseline without weather: each trajectory is the mean curve of the same
## genotype in the other seasons (held at its last value once all other
## seasons have ended), with the true initial value.
lookup_trajectory <- function(tr) {
  hold_last <- function(v) replace(v, is.nan(v), tail(v[!is.nan(v)], 1))
  P <- array(NA_real_, dim(tr$Y))
  for (y in seq_along(tr$years)) {
    m <- apply(tr$Y[, -y, , drop = FALSE], c(1, 3), mean, na.rm = TRUE)
    P[, y, ] <- t(apply(m, 1, hold_last))
    P[, y, 1] <- tr$Y[, y, 1]
  }
  P
}

## ---- Interpretation: accumulated local effects (ALE) ------------------------
## Standardised drivers on the days that occur in the data (season x day rows).
observed_points <- function(drv, tr) {
  season <- function(y) {
    n <- max(tr$nend[, y])
    matrix(drv$U[y, 1:n, ], n, dimnames = list(NULL, drv$drivers))
  }
  do.call(rbind, lapply(seq_along(drv$years), season))
}

## ---- Genotype ranking -------------------------------------------------------
## Pd, Yd: genotype x season matrices of predicted and APSIM biomass on one
## day. For each season: Spearman correlation of the two rankings and overlap
## of the predicted and true top 10 %.
rank_scores <- function(Pd, Yd, top = 0.1) {
  k <- max(1, round(top * nrow(Yd)))
  score <- function(y) {
    overlap <- intersect(order(-Pd[, y])[1:k], order(-Yd[, y])[1:k])
    data.frame(year = y, rho = cor(Pd[, y], Yd[, y], method = "spearman"),
               top = length(overlap) / k)
  }
  do.call(rbind, lapply(seq_len(ncol(Yd)), score))
}

## Baseline without weather: the genotype's mean on that day in the other
## seasons.
lookup_baseline <- function(Yd) (rowSums(Yd) - Yd) / (ncol(Yd) - 1)
