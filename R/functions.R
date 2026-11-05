## ============================================================================
##  Dynamic GxE from scratch -- all functions used in the workshop (base R only)
##
##  The pipeline is a chain of small steps. Each step takes the output of the
##  previous one, and every fitting method has the same signature, fit(A, b, ...),
##  so any step can be swapped without touching the others.
##
##    wx, traj  --make_drivers-->  drv  --make_library-->  Theta
##    traj      --make_trajectories-->  tr
##    Theta, tr --gram_stats-->  gram          (everything a linear fit needs)
##    gram      --fit_*  -->  coefficients
##    coefficients, Theta --simulate--> predicted trajectories --evaluate--> metrics
## ============================================================================

## ---- Step 1: weather -> drivers ---------------------------------------------
## Returns an array U[year, day, driver] of standardised drivers. Day 1 of the
## array is the first day in the data (das0 days after sowing).
make_drivers <- function(wx, drivers = c("TT", "DL", "RAW", "RAD")) {
  years <- sort(unique(wx$year)); das0 <- min(wx$das); nd <- max(wx$das) - das0 + 1
  U <- array(NA_real_, c(length(years), nd, length(drivers)), dimnames = list(years, NULL, drivers))
  for (i in seq_along(years)) {
    w <- wx[wx$year == years[i], ]
    U[i, w$das - das0 + 1, ] <- as.matrix(w[, drivers])
    last <- max(w$das) - das0 + 1                       # years end on different days:
    if (last < nd) U[i, (last + 1):nd, ] <- matrix(rep(U[i, last, ], each = nd - last), nd - last)  # carry last value (never used)
  }
  center <- apply(U, 3, mean); scale <- apply(U, 3, function(x) sqrt(mean((x - mean(x))^2)))
  for (j in seq_along(drivers)) U[, , j] <- (U[, , j] - center[j]) / scale[j]
  list(U = U, years = years, das0 = das0, drivers = drivers, center = center, scale = scale)
}

## ---- Step 2: drivers -> library of candidate functions -------------------------
## library_matrix(): n x d matrix of standardised drivers -> n x K matrix of candidate terms.
## Default: constant + cubic polynomial of every driver = 1 + 3 x 4 = 13 terms.
## pairs = TRUE adds the products x_i^p * x_j^q (p, q in 1:2) of every pair of drivers (37 terms).
library_matrix <- function(Z, degree = 3, pairs = FALSE) {
  pw <- function(x, p) if (p == 1) x else paste0(x, "^", p)
  d <- colnames(Z); cols <- list(rep(1, nrow(Z))); nm <- "1"
  for (j in seq_along(d)) for (p in seq_len(degree)) { cols[[length(cols) + 1]] <- Z[, j]^p; nm <- c(nm, pw(d[j], p)) }
  if (pairs) for (i in 1:(length(d) - 1)) for (j in (i + 1):length(d)) for (p in 1:2) for (q in 1:2) {
    cols[[length(cols) + 1]] <- Z[, i]^p * Z[, j]^q; nm <- c(nm, paste0(pw(d[i], p), "*", pw(d[j], q)))
  }
  L <- do.call(cbind, cols); colnames(L) <- nm; L
}
## Theta[year, day, term] = value of each candidate function on each day.
make_library <- function(drv, degree = 3, pairs = FALSE) {
  nd <- dim(drv$U); Z <- matrix(drv$U, ncol = nd[3]); colnames(Z) <- drv$drivers
  L <- library_matrix(Z, degree, pairs)
  array(L, c(nd[1:2], ncol(L)), dimnames = list(NULL, NULL, colnames(L)))
}

## ---- Step 3: long trajectory table -> array Y[genotype, year, day] -------------
## Days after maturity are NA. nend[g, y] = number of daily increments of that trajectory.
make_trajectories <- function(traj, drv, var = "biomass") {
  gen <- sort(unique(traj$genotype)); nd <- dim(drv$U)[2]
  Y <- array(NA_real_, c(length(gen), length(drv$years), nd), dimnames = list(gen, drv$years, NULL))
  Y[cbind(match(traj$genotype, gen), match(traj$year, drv$years), traj$das - drv$das0 + 1)] <- traj[[var]]
  nend <- apply(Y, c(1, 2), function(v) sum(!is.na(v)) - 1)
  list(Y = Y, nend = nend, genotypes = gen, years = drv$years, das0 = drv$das0)
}

## ---- Step 4: the regression problem ------------------------------------------
## For each trajectory the daily growth increments are  z_t = y_{t+1} - y_t  and the
## regressors are the library values on day t:  X_t = Theta[year, t, ].  Stacking all
## days of all trajectories gives the design matrix; build_design() shows it explicitly.
build_design <- function(Theta, tr, genotypes = tr$genotypes) {
  X <- list(); z <- list(); id <- list()
  for (g in genotypes) for (y in seq_along(tr$years)) {
    gi <- match(g, tr$genotypes); n <- tr$nend[gi, y]; v <- tr$Y[gi, y, 1:(n + 1)]
    X[[length(X) + 1]] <- matrix(Theta[y, 1:n, ], n); z[[length(z) + 1]] <- diff(v)
    id[[length(id) + 1]] <- data.frame(genotype = g, year = tr$years[y], day = 0:(n - 1))
  }
  X <- do.call(rbind, X); colnames(X) <- dimnames(Theta)[[3]]
  list(X = X, z = unlist(z), id = do.call(rbind, id))
}

## Everything a linear regression needs from the data is  A = X'X  and  b = X'z.
## We store them per genotype and per year, so that leave-one-year-out cross-validation
## only needs additions and subtractions of these small matrices (no refitting from scratch).
gram_stats <- function(Theta, tr) {
  G <- dim(tr$Y)[1]; ny <- length(tr$years); K <- dim(Theta)[3]
  A <- array(0, c(G, ny, K, K)); b <- array(0, c(G, ny, K))
  for (g in seq_len(G)) for (y in seq_len(ny)) {
    n <- tr$nend[g, y]; X <- matrix(Theta[y, 1:n, ], n)
    A[g, y, , ] <- crossprod(X); b[g, y, ] <- crossprod(X, diff(tr$Y[g, y, 1:(n + 1)]))
  }
  list(A = A, b = b, n = tr$nend, names = dimnames(Theta)[[3]])
}

## Sum the statistics over a set of years (the training years of one fold).
sum_years <- function(gram, years) {
  list(A = apply(gram$A[, years, , , drop = FALSE], c(1, 3, 4), sum),    # genotype x K x K
       b = apply(gram$b[, years, , drop = FALSE], c(1, 3), sum),          # genotype x K
       n = apply(gram$n[, years, drop = FALSE], 1, sum))                  # genotype
}

## Put all terms on a common scale (unit root-mean-square column), so that coefficients
## can be compared and penalised fairly. Returns scaled A, b and the scale s.
scale_gram <- function(A, b, n, s = NULL) {
  if (is.null(s)) s <- sqrt(diag(A) / n)                       # RMS of every library column
  list(A = A / n / outer(s, s), b = b / n / s, s = s)
}

## ---- Step 5: fitting methods (all with the same signature fit(A, b, ...)) ----------
## Every method returns list(coef, support) with coef on the scaled problem.
fit_ols <- function(A, b, ...) {
  list(coef = solve(A, b), support = seq_along(b))
}
fit_ridge <- function(A, b, lambda = 0.01, ...) {
  list(coef = solve(A + lambda * mean(diag(A)) * diag(nrow(A)), b), support = seq_along(b))
}
## Sequential thresholding / backward elimination: refit, drop the smallest coefficient,
## repeat until k terms remain.
fit_ssr <- function(A, b, k = 8, alpha = 1e-6, ...) {
  S <- seq_along(b)
  repeat {
    co <- solve(A[S, S, drop = FALSE] + alpha * diag(length(S)), b[S])
    if (length(S) <= k) break
    S <- S[-which.min(abs(co))]
  }
  coef <- numeric(length(b)); coef[S] <- co
  list(coef = coef, support = S)
}

## ---- Step 5b: shared coefficients + genotype deviations -----------------------------
## xi_g = mu + delta_g on a fixed support: mu from the pooled problem (light ridge),
## delta_g from a ridge-shrunk fit of each genotype's own residual growth.
fit_hier <- function(Ag, bg, support, lam_mu = 1e-4, lam_g = 0.01) {
  K <- dim(Ag)[2]; k <- length(support); G <- dim(Ag)[1]
  Apool <- apply(Ag[, support, support, drop = FALSE], c(2, 3), sum); bpool <- colSums(bg[, support, drop = FALSE])
  mu <- numeric(K); mu[support] <- solve(Apool + lam_mu * mean(diag(Apool)) * diag(k), bpool)
  delta <- matrix(0, G, K)
  for (g in seq_len(G)) {
    Aq <- Ag[g, support, support]; tr <- sum(diag(Aq)) / k
    delta[g, support] <- solve(Aq + lam_g * tr * diag(k), bg[g, support] - Aq %*% mu[support])
  }
  list(mu = mu, delta = delta, xi = sweep(delta, 2, mu, "+"), support = support)
}

## ---- Step 6: simulate (roll the equation forward) ----------------------------------
## For a library that only depends on the weather the "ODE solution" is a running total:
##   y(t) = y0 + sum_{s < t} Theta(s) %*% xi        -> no ODE solver needed.
## xi: genotype x term matrix (unscaled coefficients); Theta_y: day x term for one year.
simulate <- function(xi, Theta_y, y0) {
  inc <- Theta_y[-nrow(Theta_y), , drop = FALSE] %*% t(xi)         # (days-1) x genotype increments
  sweep(rbind(0, apply(inc, 2, cumsum)), 2, y0, "+")               # day x genotype
}

## ---- Step 7: leave-one-year-out cross-validation -------------------------------------
## fitter: any function fit(A, b, ...) returning list(coef, support); extra arguments are passed on.
## support_k: number of library terms kept (used to pick the shared support by `fitter`).
## Fit on a set of years: pooled fit by `fitter` chooses the support, then shared + genotype deviations.
## Returns coefficients in the original (unscaled) units of the library.
fit_years <- function(gram, years, fitter = fit_ssr, k = 12, lam_mu = 1e-4, lam_g = 0.01, ...) {
  K <- dim(gram$A)[3]; st <- sum_years(gram, years); ntot <- sum(st$n)
  Apool <- apply(st$A, c(2, 3), sum); bpool <- colSums(st$b); s <- sqrt(diag(Apool) / ntot)
  Ag <- st$A / ntot; for (i in 1:K) for (j in 1:K) Ag[, i, j] <- st$A[, i, j] / ntot / (s[i] * s[j])
  bg <- sweep(st$b / ntot, 2, s, "/")
  sp <- fitter(Apool / ntot / outer(s, s), bpool / ntot / s, k = k, ...)$support
  h <- fit_hier(Ag, bg, sp, lam_mu, lam_g)
  list(xi = sweep(h$xi, 2, s, "/"), mu = h$mu / s, support = sp, names = gram$names)
}

cv_blocked <- function(Theta, tr, gram, fitter = fit_ssr, k = 12, lam_mu = 1e-4, lam_g = 0.01, ...) {
  ny <- length(tr$years); G <- length(tr$genotypes); nd <- dim(Theta)[2]
  P <- array(NA_real_, c(G, ny, nd)); S_all <- vector("list", ny)
  for (yv in seq_len(ny)) {
    train <- setdiff(seq_len(ny), yv); f <- fit_years(gram, train, fitter, k, lam_mu, lam_g, ...)
    y0 <- apply(tr$Y[, train, 1, drop = FALSE], 1, mean)                        # initial biomass: mean of training years
    P[, yv, ] <- t(simulate(f$xi, Theta[yv, , ], y0)); S_all[[yv]] <- f$support
  }
  list(pred = P, support = S_all)
}

## ---- Step 8: metrics ------------------------------------------------------------------
## RMSE over the whole trajectory (days up to maturity), at a fixed day, and at maturity.
evaluate <- function(P, tr, fixed_day = 99) {
  G <- dim(P)[1]; ny <- dim(P)[2]; se <- 0; cnt <- 0; fin <- numeric(0)
  for (g in seq_len(G)) for (y in seq_len(ny)) {
    n <- tr$nend[g, y]; d <- 2:(n + 1)                      # day index 1..n (day 0 is the given start)
    se <- se + sum((P[g, y, d] - tr$Y[g, y, d])^2); cnt <- cnt + length(d)
    fin <- c(fin, P[g, y, n + 1] - tr$Y[g, y, n + 1])
  }
  c(trajectory_rmse = sqrt(se / cnt), fixed_day_rmse = sqrt(mean((P[, , fixed_day] - tr$Y[, , fixed_day])^2, na.rm = TRUE)),
    maturity_rmse = sqrt(mean(fin^2)))
}

## ---- Step 9: reading the model --------------------------------------------------------
## Standardised drivers on the days that actually occur in the data (rows = year x day).
observed_points <- function(drv, tr) {
  Z <- NULL
  for (y in seq_along(drv$years)) {
    n <- max(tr$nend[, y]); Z <- rbind(Z, matrix(drv$U[y, 1:n, ], n, dimnames = list(NULL, drv$drivers)))
  }
  Z
}
## Accumulated local effects (ALE) of driver j on the growth rate. We move only driver j across a
## small interval around each observed value and average the change in the rate. Unlike partial
## dependence it never asks the model about weather combinations that do not occur.
## rate: function(Z) -> growth rate for an n x d matrix of standardised drivers.
ale_curve <- function(rate, Z, j, K = 24) {
  e <- unique(quantile(Z[, j], seq(0, 1, length.out = K + 1))); nb_ <- length(e) - 1
  b <- pmin(pmax(findInterval(Z[, j], e), 1), nb_)
  Zu <- Z; Zl <- Z; Zu[, j] <- e[b + 1]; Zl[, j] <- e[b]
  d <- rate(Zu) - rate(Zl); nb <- tabulate(b, nb_)
  md <- vapply(split(d, factor(b, levels = seq_len(nb_))), sum, 0) / pmax(nb, 1)
  A <- c(0, cumsum(md)); mid <- (A[-1] + A[-length(A)]) / 2
  data.frame(z = e, effect = A - sum(mid * nb) / sum(nb))
}

## ---- Step 10: ranking genotypes on a fixed day ----------------------------------------------
## Pd, Yd: genotype x year matrices (predicted and APSIM biomass on the same day).
## For every year: Spearman correlation between predicted and true ranking of the genotypes,
## and the overlap of the predicted and true top 10 %.
rank_scores <- function(Pd, Yd, top = 0.1) {
  k <- max(1, round(top * nrow(Yd)))
  do.call(rbind, lapply(seq_len(ncol(Yd)), function(y) data.frame(
    year = y, rho = cor(Pd[, y], Yd[, y], method = "spearman"),
    top = length(intersect(order(-Pd[, y])[1:k], order(-Yd[, y])[1:k])) / k)))
}
## Baseline that uses no weather at all: the genotype's mean on that day in the other years.
lookup_baseline <- function(Yd) (rowSums(Yd) - Yd) / (ncol(Yd) - 1)
