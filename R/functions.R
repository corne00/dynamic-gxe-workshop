## =============================================================================
##  Sparse identification of crop growth: functions used in the workshop
##
##    weather --make_drivers--> drv --make_library--> Theta   (candidate terms)
##    biomass --make_trajectories--> tr
##    Theta[year, day, term], tr --make_design--> Theta (rows = days), dy
##    fit_*(Theta, dy) --> coefficients --simulate--> prediction --evaluate--> RMSE
##
##  Base R, except for glmnet and leaps (fitting).
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
## 1..degree of every driver; pairs = TRUE adds products of pairs of drivers
## with exponents p, q >= 1 and p + q <= degree.
## The constant term is the intercept of the regression and is not included.
library_matrix <- function(Z, degree = 3, pairs = FALSE) {
  term <- function(x, p) if (p == 1) x else paste0(x, "^", p)
  cols <- list()
  for (d in colnames(Z)) {
    for (p in 1:degree) cols[[term(d, p)]] <- Z[, d]^p
  }
  if (pairs) {
    for (ij in combn(colnames(Z), 2, simplify = FALSE)) {
      for (p in seq_len(degree - 1)) {
        for (q in seq_len(degree - p)) {
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
## Rows: days 2, ..., n + 1 of every year (all days with a record before them).
## Biomass is recorded at the end of a day, so day i drives the change from
## observation i - 1 to i. Row i of the returned matrix Theta is the library of
## day i, and dy[i] = y_i - y_(i-1) is the growth rate at the middle of day i.
make_design <- function(Theta, tr, genotype) {
  g <- match(genotype, tr$genotypes)
  n <- tr$nend[g, ]
  years <- seq_along(tr$years)
  rows <- lapply(years, function(y) Theta[y, 2:(n[y] + 1), ])
  dy <- lapply(years, function(y) diff(tr$Y[g, y, 1:(n[y] + 1)]))
  rows <- do.call(rbind, rows)
  colnames(rows) <- dimnames(Theta)[[3]]
  list(Theta = rows, dy = unlist(dy), year = rep(years, n))
}

## ---- Fitting methods: fit(Theta, dy) -> coefficients (intercept first) ----------
named <- function(b, Theta) setNames(unname(b), c("1", colnames(Theta)))

col_sd <- function(Theta) {
  n <- nrow(Theta)
  sqrt((colSums(Theta^2) - n * colMeans(Theta)^2) / (n - 1))
}

fit_ols <- function(Theta, dy) {
  Z <- cbind(`1` = 1, Theta)
  drop(solve(crossprod(Z), crossprod(Z, dy)))
}

## Ridge in closed form, with the penalty scaling of glmnet.
fit_ridge <- function(Theta, dy, lambda = 1e-3) {
  n <- nrow(Theta); mx <- colMeans(Theta)
  var_x <- col_sd(Theta)^2 * (n - 1) / n
  sy <- sd(dy) * sqrt((n - 1) / n)
  b <- solve(crossprod(Theta) / n - tcrossprod(mx) + diag(lambda / sy * var_x),
             crossprod(Theta, dy) / n - mx * mean(dy))
  c(`1` = mean(dy) - sum(mx * b), setNames(drop(b), colnames(Theta)))
}

## Lasso and elastic net by glmnet.
fit_glmnet <- function(Theta, dy, alpha, lambda) {
  fit <- glmnet::glmnet(Theta, dy, alpha = alpha, lambda = lambda,
                        type.gaussian = "covariance")
  named(as.vector(coef(fit)), Theta)
}
fit_lasso <- function(Theta, dy, lambda = 3e-4) fit_glmnet(Theta, dy, 1, lambda)
fit_enet <- function(Theta, dy, lambda = 6e-4) fit_glmnet(Theta, dy, 0.5, lambda)

## Best subset or backward elimination with k terms.
fit_subset <- function(Theta, dy, k = 8, method = "exhaustive") {
  b <- coef(leaps::regsubsets(Theta, dy, nvmax = k, method = method), k)
  full <- setNames(numeric(ncol(Theta) + 1), c("1", colnames(Theta)))
  full[c("1", names(b)[-1])] <- b
  full
}

## Sequential thresholded least squares: refit, drop terms below the threshold.
fit_stlsq <- function(Theta, dy, threshold = 0.05) {
  Z <- cbind(`1` = 1, Theta)
  G <- crossprod(Z)
  rhs <- drop(crossprod(Z, dy))
  size <- col_sd(Theta) / sd(dy)
  active <- colnames(Theta)
  repeat {
    i <- c("1", active)
    b <- drop(solve(G[i, i], rhs[i]))
    keep <- abs(b[-1]) * size[active] >= threshold
    if (all(keep) || !any(keep)) break
    active <- active[keep]
  }
  full <- setNames(numeric(ncol(Theta) + 1), c("1", colnames(Theta)))
  full[names(b)] <- b
  full
}

## Methods compared in the workshop (hyperparameters fixed for simplicity).
fit_methods <- list(
  OLS = fit_ols,
  ridge = fit_ridge,
  lasso = fit_lasso,
  `elastic net` = fit_enet,
  `best subset (k = 8)` = function(Theta, dy) fit_subset(Theta, dy, k = 8),
  `backward (k = 8)` = function(Theta, dy) fit_subset(Theta, dy, 8, "backward"),
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
      fit(d$Theta[train, ], d$dy[train])
    }))
    y0 <- (sum(tr$Y[g, , 1]) - tr$Y[g, , 1]) / (S - 1)
    P[g, , ] <- t(simulate(coef, Theta, y0))
  }
  P
}

## ---- State feedback and coupled models --------------------------------------
## The modelled states enter the library through standardised terms (centre and
## standard deviation of the observed values of the genotype): powers of each
## state, and with products = TRUE also products of states and products of a
## state with the powers 1 and 2 of each driver. trs is a list of trajectory
## sets from make_trajectories(), one per modelled state, named by state; every
## state has its own equation, and all are fitted with the same fit(Theta, dy). The
## state terms of a row use the mean of the records on either side of the day.
state_scale <- function(Y) {
  m <- mean(Y, na.rm = TRUE)
  c(center = m, scale = sqrt(mean((Y - m)^2, na.rm = TRUE)))
}

## Library columns that interact with a state: drivers and their squares.
interacting_terms <- function(Theta, drivers) {
  want <- unlist(lapply(drivers, function(d) c(d, paste0(d, "^2"))))
  want[want %in% dimnames(Theta)[[3]]]
}

## z: rows x states (standardised, named), U: rows x interacting library columns.
## Products have total exponent at most degree.
state_terms <- function(z, U, degree = 3, pairs = FALSE) {
  term <- function(x, p) if (p == 1) x else paste0(x, "^", p)
  cols <- list(); st <- colnames(z)
  for (a in st) for (p in seq_len(degree)) cols[[term(a, p)]] <- z[, a]^p
  if (pairs) {
    if (length(st) > 1) for (ij in combn(st, 2, simplify = FALSE)) {
      for (p in seq_len(degree - 1)) for (q in seq_len(degree - p)) {
        cols[[paste0(term(ij[1], p), "*", term(ij[2], q))]] <- z[, ij[1]]^p * z[, ij[2]]^q
      }
    }
    power <- function(u) if (grepl("^", u, fixed = TRUE)) as.integer(sub(".*\\^", "", u)) else 1L
    for (a in st) for (u in colnames(U)) for (p in seq_len(degree - 1)) {
      if (p + power(u) <= degree) cols[[paste0(term(a, p), "*", u)]] <- z[, a]^p * U[, u]
    }
  }
  do.call(cbind, cols)
}

## Regression problem of one genotype with state terms: one column of dy per
## state; the state terms use the state at the observation of each row.
make_design_states <- function(Theta, trs, genotype, drivers, degree = 3,
                               pairs = FALSE, feedback = TRUE) {
  d <- make_design(Theta, trs[[1]], genotype)
  g <- match(genotype, trs[[1]]$genotypes); n <- trs[[1]]$nend[g, ]
  dy <- sapply(trs, function(tr) make_design(Theta, tr, genotype)$dy)
  scales <- sapply(trs, function(tr) state_scale(tr$Y[g, , ]))
  full <- d$Theta
  if (feedback) {
    z <- sapply(names(trs), function(k) {
      y <- unlist(lapply(seq_along(n), function(i) {   # state in the middle of the day
        Y <- trs[[k]]$Y[g, i, 1:(n[i] + 1)]
        (Y[-length(Y)] + Y[-1]) / 2
      }))
      (y - scales["center", k]) / scales["scale", k]
    })
    U <- d$Theta[, interacting_terms(Theta, drivers), drop = FALSE]
    full <- cbind(full, state_terms(matrix(z, ncol = length(trs), dimnames = list(NULL, names(trs))),
                                    U, degree, pairs))
  }
  list(Theta = full, dy = matrix(dy, ncol = length(trs), dimnames = list(NULL, names(trs))),
       year = d$year, scales = scales)
}

## Solves the coupled system dy_k/dt = (1, theta(y, u)) coef_k. coef: list of
## year x (1 + terms) matrices, one per state; y0: year x state. Returns a list
## of day x year matrices. terms(y, u) is the library, with y the year x state
## matrix. State terms can make a solution blow up; the years that do are solved
## one by one and returned as NA (a state beyond limit[k] counts as blown up), so
## that the others are not lost.
simulate_states <- function(coef, U, y0, terms = function(y, u) u, method = "lsoda",
                            limit = rep(1e4, length(coef))) {
  D <- dim(U)[2]; m <- length(coef)
  solve_years <- function(idx) {
    S <- length(idx)
    days <- lapply(seq_len(D), function(j) matrix(U[idx, j, ], S))
    B <- lapply(coef, function(b) b[idx, , drop = FALSE])
    rhs <- function(t, y, parms) {
      pos <- min(max(t + 0.5, 0), D - 1)
      j <- min(floor(pos), D - 2)
      w <- pos - j
      u <- (1 - w) * days[[j + 1]] + w * days[[j + 2]]
      Tm <- terms(matrix(y, S, m), u)
      list(unlist(lapply(seq_len(m), function(k) {
        B[[k]][, 1] + .rowSums(Tm * B[[k]][, -1, drop = FALSE], S, ncol(Tm))
      })))
    }
    out <- tryCatch(suppressWarnings(deSolve::ode(
      as.vector(y0[idx, , drop = FALSE]), 0:(D - 1), rhs, NULL, method = method,
      rtol = 1e-4, atol = 1e-4)[, -1, drop = FALSE]), error = function(e) NULL)
    if (is.null(out) || nrow(out) < D || any(!is.finite(out))) return(NULL)
    for (k in seq_len(m)) if (any(abs(out[, (k - 1) * S + seq_len(S)]) > limit[k])) return(NULL)
    out
  }
  S <- dim(U)[1]
  pack <- function(out, idx) lapply(seq_len(m), function(k) {
    out[, (k - 1) * length(idx) + seq_along(idx), drop = FALSE]
  })
  all <- solve_years(seq_len(S))
  if (!is.null(all)) return(pack(all, seq_len(S)))
  res <- lapply(seq_len(m), function(k) matrix(NA_real_, D, S))
  for (s in seq_len(S)) {
    one <- solve_years(s)
    if (is.null(one)) next
    for (k in seq_len(m)) res[[k]][, s] <- one[, k]
  }
  res
}

## Leave-one-year-out validation of the (coupled) model: a list of P[genotype,
## year, day], one per state.
cross_validate_states <- function(Theta, trs, fit, genotypes = trs[[1]]$genotypes,
                                  drivers, degree = 3, pairs = FALSE,
                                  feedback = TRUE) {
  tr <- trs[[1]]; S <- length(tr$years); m <- length(trs)
  P <- lapply(trs, function(t) array(NA_real_, dim(t$Y))); diverged <- 0
  inter <- interacting_terms(Theta, drivers)
  inter_at <- match(inter, dimnames(Theta)[[3]])
  for (g in match(genotypes, tr$genotypes)) {
    d <- make_design_states(Theta, trs, tr$genotypes[g], drivers, degree, pairs, feedback)
    coef <- lapply(seq_len(m), function(k) t(sapply(seq_len(S), function(v) {
      train <- d$year != v
      fit(d$Theta[train, ], d$dy[train, k])
    })))
    y0 <- sapply(trs, function(t) (sum(t$Y[g, , 1]) - t$Y[g, , 1]) / (S - 1))
    terms <- function(y, u) u
    if (feedback) terms <- function(y, u) {
      z <- sweep(sweep(y, 2, d$scales["center", ]), 2, d$scales["scale", ], "/")
      colnames(z) <- names(trs)
      ui <- u[, inter_at, drop = FALSE]; colnames(ui) <- inter
      cbind(u, state_terms(z, ui, degree, pairs))
    }
    limit <- 10 * sapply(trs, function(t) max(abs(t$Y[g, , ]), na.rm = TRUE))
    sim <- simulate_states(coef, Theta, matrix(y0, S, m), terms, limit = limit)
    for (k in seq_len(m)) P[[k]][g, , ] <- t(sim[[k]])
    diverged <- diverged + sum(is.na(sim[[1]][1, ]))
  }
  attr(P, "diverged") <- diverged   # held-out years whose solution blew up (NA)
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

