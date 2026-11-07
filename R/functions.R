## =============================================================================
##  Sparse identification of crop growth: functions used in the workshop
##
##    weather --make_drivers--> drv --make_library--> Theta   (candidate terms)
##    biomass --make_trajectories--> tr
##    Theta, tr --make_design--> X, dy                (growth rate by centered differences)
##    fit_*(X, dy) --> coefficients --simulate--> prediction --evaluate--> RMSE
##
##  Base R, except for pracma (derivatives), and glmnet and leaps (fitting).
## =============================================================================

## ---- Drivers: array U[year, day, driver], standardised -------------------
make_drivers <- function(wx, drivers = c("Tmin", "Tmax", "Rain", "RAD", "VP", "DL",
                                          "TT", "RAW")) {
  years <- sort(unique(wx$year))
  das0 <- min(wx$das)
  nd <- max(wx$das) - das0 + 1
  U <- array(NA_real_, c(length(years), nd, length(drivers)),
             dimnames = list(years, NULL, drivers))
  for (j in seq_along(drivers)) {
    i <- cbind(match(wx$year, years), wx$das - das0 + 1, j)
    U[i] <- wx[[drivers[j]]]
  }
  # pad years that end early with their last value
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

## Theta[year, day, term]
make_library <- function(drv, degree = 3, pairs = FALSE) {
  d <- dim(drv$U)
  Z <- matrix(drv$U, ncol = d[3], dimnames = list(NULL, drv$drivers))
  L <- library_matrix(Z, degree, pairs)
  array(L, c(d[1:2], ncol(L)), dimnames = list(NULL, NULL, colnames(L)))
}

## ---- Extra terms: functions of the raw weather, standardised like the drivers
## exprs: named character vector in R syntax. Available are the columns of wx,
## Tmean, VPD (hPa) and das, and the functions roll(x, k) (mean of the last k
## days), cum(x) and lag(x, k); max and min act elementwise.
extra_terms <- function(wx, exprs, years = sort(unique(wx$year))) {
  das0 <- min(wx$das); nd <- max(wx$das) - das0 + 1
  fun <- new.env(parent = baseenv())
  fun$max <- pmax; fun$min <- pmin; fun$cum <- cumsum
  fun$roll <- function(x, k) {
    vapply(seq_along(x), function(t) mean(x[max(1, t - k + 1):t]), 0)
  }
  fun$lag <- function(x, k) x[pmax(seq_along(x) - k, 1)]
  X <- array(NA_real_, c(length(years), nd, length(exprs)),
             dimnames = list(years, NULL, names(exprs)))
  for (y in seq_along(years)) {
    w <- wx[wx$year == years[y], ]
    v <- lapply(w[setdiff(names(w), c("year", "das"))], function(x) {
      x <- c(x, rep(NA, nd - length(x)))
      x[is.na(x)] <- tail(x[!is.na(x)], 1)     # later days repeat the last value
      x
    })
    v$das <- das0 + seq_len(nd) - 1
    v$Tmean <- (v$Tmin + v$Tmax) / 2
    v$VPD <- 6.108 * exp(17.27 * v$Tmax / (v$Tmax + 237.3)) - v$VP
    for (j in seq_along(exprs)) X[y, , j] <- eval(parse(text = exprs[[j]]), v, fun)
  }
  for (j in seq_along(exprs)) {
    x <- X[, , j]
    X[, , j] <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  }
  X
}

## Appends extra terms to a library Theta[year, day, term].
bind_terms <- function(Theta, extra) {
  d <- dim(Theta)
  array(c(Theta, extra), c(d[1:2], d[3] + dim(extra)[3]),
        dimnames = list(NULL, NULL, c(dimnames(Theta)[[3]], dimnames(extra)[[3]])))
}

## ---- Measured traits as drivers ---------------------------------------------
## Adds measured states of one genotype (columns of traj, for example "lai" or
## "biomass") to the weather table, so that they can be used like drivers. A
## record is the value at the end of a day, so the value of a day is the mean of
## the records at its start and end.
with_traits <- function(wx, traj, genotype, vars = "lai") {
  t <- traj[traj$genotype == genotype, ]
  for (v in vars) {
    wx[[v]] <- NA_real_
    for (y in unique(wx$year)) {
      r <- t[t$year == y, ]
      r <- r[order(r$das), ]
      day <- c(r[[v]][1], (r[[v]][-1] + r[[v]][-nrow(r)]) / 2)
      i <- wx$year == y
      wx[[v]][i] <- day[match(wx$das[i], r$das)]
    }
  }
  wx
}

## ---- Trajectories: array Y[genotype, year, day], NA after maturity --------
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
## Rows: observations 2, ..., n of every year. dy: centered difference of the
## biomass. Biomass is recorded at the end of a day, so day i drives the change
## from observation i - 1 to i. X: library at observation i, the mean of the
## library of days i and i + 1.
make_design <- function(Theta, tr, genotype) {
  g <- match(genotype, tr$genotypes)
  n <- tr$nend[g, ]
  years <- seq_along(tr$years)
  X <- lapply(years, function(y)
    (Theta[y, 2:n[y], ] + Theta[y, 3:(n[y] + 1), ]) / 2)
  dy <- lapply(years, function(y)
    pracma::gradient(tr$Y[g, y, 1:(n[y] + 1)])[2:n[y]])
  X <- do.call(rbind, X)
  colnames(X) <- dimnames(Theta)[[3]]
  list(X = X, dy = unlist(dy), year = rep(years, n - 1))
}

## ---- Fitting methods: fit(X, dy) -> coefficients (intercept first) ----------
named <- function(b, X) setNames(unname(b), c("1", colnames(X)))

col_sd <- function(X) {
  n <- nrow(X)
  sqrt((colSums(X^2) - n * colMeans(X)^2) / (n - 1))
}

fit_ols <- function(X, dy) {
  Z <- cbind(`1` = 1, X)
  drop(solve(crossprod(Z), crossprod(Z, dy)))
}

## Ridge in closed form, with the penalty scaling of glmnet.
fit_ridge <- function(X, dy, lambda = 1e-3) {
  n <- nrow(X); mx <- colMeans(X)
  var_x <- col_sd(X)^2 * (n - 1) / n
  sy <- sd(dy) * sqrt((n - 1) / n)
  b <- solve(crossprod(X) / n - tcrossprod(mx) + diag(lambda / sy * var_x),
             crossprod(X, dy) / n - mx * mean(dy))
  c(`1` = mean(dy) - sum(mx * b), setNames(drop(b), colnames(X)))
}

## Lasso and elastic net by glmnet.
fit_glmnet <- function(X, dy, alpha, lambda) {
  fit <- glmnet::glmnet(X, dy, alpha = alpha, lambda = lambda,
                        type.gaussian = "covariance")
  named(as.vector(coef(fit)), X)
}
fit_lasso <- function(X, dy, lambda = 3e-4) fit_glmnet(X, dy, 1, lambda)
fit_enet <- function(X, dy, lambda = 6e-4) fit_glmnet(X, dy, 0.5, lambda)

## Best subset or backward elimination with k terms.
fit_subset <- function(X, dy, k = 8, method = "exhaustive") {
  b <- coef(leaps::regsubsets(X, dy, nvmax = k, method = method), k)
  full <- setNames(numeric(ncol(X) + 1), c("1", colnames(X)))
  full[c("1", names(b)[-1])] <- b
  full
}

## Sequential thresholded least squares: refit, drop terms below the threshold.
fit_stlsq <- function(X, dy, threshold = 0.05) {
  Z <- cbind(`1` = 1, X)
  G <- crossprod(Z)
  rhs <- drop(crossprod(Z, dy))
  size <- col_sd(X) / sd(dy)
  active <- colnames(X)
  repeat {
    i <- c("1", active)
    b <- drop(solve(G[i, i], rhs[i]))
    keep <- abs(b[-1]) * size[active] >= threshold
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

## ---- Simulation -------------------------------------------------------------
## Solves dy/dt = (1, theta(y, u)) coef. coef: year x (1 + terms), U: year
## x day x variable. terms(y, u) is the library; by default U is the library.
simulate <- function(coef, U, y0, terms = function(y, u) u, method = "lsoda") {
  S <- dim(U)[1]; D <- dim(U)[2]
  days <- lapply(seq_len(D), function(j) matrix(U[, j, ], S))
  b0 <- coef[, 1]; B <- coef[, -1, drop = FALSE]
  rhs <- function(t, y, parms) {
    pos <- min(max(t + 0.5, 0), D - 1)      # day j is placed at time j - 3/2
    j <- min(floor(pos), D - 2)
    w <- pos - j                            # linear interpolation between days
    u <- (1 - w) * days[[j + 1]] + w * days[[j + 2]]
    list(b0 + .rowSums(terms(y, u) * B, S, ncol(B)))
  }
  deSolve::ode(y0, 0:(D - 1), rhs, NULL, method = method, rtol = 1e-4,
               atol = 1e-4)[, -1, drop = FALSE]
}

## n genotypes spread evenly over the ranking of mean final biomass.
spread_genotypes <- function(tr, n = 3) {
  last <- cbind(as.vector(row(tr$nend)), as.vector(col(tr$nend)),
                as.vector(tr$nend) + 1)
  final <- matrix(tr$Y[last], nrow(tr$nend), dimnames = dimnames(tr$nend))
  ranked <- names(sort(rowMeans(final)))
  ranked[round(seq(1, length(ranked), length.out = n))]
}

## ---- Validation: leave-one-year-out ---------------------------------------
## Each year is predicted by a model fitted to the other years of the same
## genotype, starting from the mean initial biomass of those years.
## Returns P[genotype, year, day].
cross_validate <- function(Theta, tr, fit, genotypes = tr$genotypes) {
  P <- array(NA_real_, dim(tr$Y))
  S <- length(tr$years)
  for (g in match(genotypes, tr$genotypes)) {
    d <- make_design(Theta, tr, tr$genotypes[g])
    coef <- t(sapply(seq_len(S), function(v) {
      train <- d$year != v
      fit(d$X[train, ], d$dy[train])
    }))
    y0 <- (sum(tr$Y[g, , 1]) - tr$Y[g, , 1]) / (S - 1)
    P[g, , ] <- t(simulate(coef, Theta, y0))
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

## ---- Interpretation: accumulated local effects (ALE) ------------------------
## Standardised drivers on the days that occur in the data (year x day rows).
observed_points <- function(drv, tr) {
  year <- function(y) {
    n <- max(tr$nend[, y])
    matrix(drv$U[y, 1:n, ], n, dimnames = list(NULL, drv$drivers))
  }
  do.call(rbind, lapply(seq_along(drv$years), year))
}

## ---- Genotype ranking -------------------------------------------------------
## Pd, Yd: genotype x year matrices of predicted and APSIM biomass on one
## day. For each year: Spearman correlation of the two rankings and overlap
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

