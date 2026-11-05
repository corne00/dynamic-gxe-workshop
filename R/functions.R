## ============================================================================
##  Sparse identification of crop growth: all functions used in the workshop (base R)
##
##    wx   --make_drivers-->       drv --make_library--> Theta     (candidate functions)
##    traj --make_trajectories-->  tr                              (biomass trajectories)
##    Theta, tr --gram_stats-->    gram                            (A = X'X and b = X'z)
##    gram --fit_*-->              coefficients --simulate--> predictions --evaluate--> RMSE
##
##  Every fitting method has the signature fit(A, b, ...), so methods are interchangeable.
## ============================================================================

## ---- Drivers: array U[year, day, driver], standardised ---------------------------
make_drivers <- function(wx, drivers = c("TT", "DL", "RAW", "RAD")) {
  years <- sort(unique(wx$year)); das0 <- min(wx$das); nd <- max(wx$das) - das0 + 1
  U <- array(NA_real_, c(length(years), nd, length(drivers)), dimnames = list(years, NULL, drivers))
  for (j in seq_along(drivers)) U[cbind(match(wx$year, years), wx$das - das0 + 1, j)] <- wx[[drivers[j]]]
  for (y in seq_along(years)) for (j in seq_along(drivers)) {      # seasons end on different days:
    u <- U[y, , j]; U[y, is.na(u), j] <- tail(u[!is.na(u)], 1)     # pad with the last value (never used in fits)
  }
  center <- apply(U, 3, mean)
  scale  <- apply(U, 3, function(x) sqrt(mean((x - mean(x))^2)))
  list(U = sweep(sweep(U, 3, center), 3, scale, "/"), years = years, das0 = das0,
       drivers = drivers, center = center, scale = scale)
}

## ---- Library of candidate functions ------------------------------------------------
## Z: n x d matrix of standardised drivers -> n x K matrix of terms.
## A constant and powers 1..degree of every driver; pairs = TRUE adds products of pairs.
library_matrix <- function(Z, degree = 3, pairs = FALSE) {
  term <- function(x, p) if (p == 1) x else paste0(x, "^", p)
  cols <- list(`1` = rep(1, nrow(Z)))
  for (d in colnames(Z)) for (p in 1:degree) cols[[term(d, p)]] <- Z[, d]^p
  if (pairs) for (ij in combn(colnames(Z), 2, simplify = FALSE)) for (p in 1:2) for (q in 1:2)
    cols[[paste0(term(ij[1], p), "*", term(ij[2], q))]] <- Z[, ij[1]]^p * Z[, ij[2]]^q
  do.call(cbind, cols)
}
## Theta[year, day, term]
make_library <- function(drv, degree = 3, pairs = FALSE) {
  d <- dim(drv$U); Z <- matrix(drv$U, ncol = d[3], dimnames = list(NULL, drv$drivers))
  L <- library_matrix(Z, degree, pairs)
  array(L, c(d[1:2], ncol(L)), dimnames = list(NULL, NULL, colnames(L)))
}

## ---- Trajectories: array Y[genotype, year, day], NA after maturity -----------------
## nend[g, y] is the number of daily steps of trajectory (g, y).
make_trajectories <- function(traj, drv, var = "biomass") {
  gen <- sort(unique(traj$genotype)); nd <- dim(drv$U)[2]
  Y <- array(NA_real_, c(length(gen), length(drv$years), nd), dimnames = list(gen, drv$years, NULL))
  Y[cbind(match(traj$genotype, gen), match(traj$year, drv$years), traj$das - drv$das0 + 1)] <- traj[[var]]
  nend <- apply(Y, c(1, 2), function(v) sum(!is.na(v)) - 1)
  list(Y = Y, nend = nend, genotypes = gen, years = drv$years, das0 = drv$das0)
}

## ---- Regression problem ---------------------------------------------------------------
## One trajectory: increments z_t = y_{t+1} - y_t regressed on the library at day t.
trajectory_xz <- function(Theta, tr, g, y) {
  n <- tr$nend[g, y]
  list(X = matrix(Theta[y, 1:n, ], n), z = diff(tr$Y[g, y, 1:(n + 1)]))
}

## Explicit design matrix X and response z for a set of genotypes (all seasons stacked).
build_design <- function(Theta, tr, genotypes = tr$genotypes) {
  idx <- expand.grid(y = seq_along(tr$years), g = match(genotypes, tr$genotypes))
  xz  <- Map(function(g, y) trajectory_xz(Theta, tr, g, y), idx$g, idx$y)
  n   <- tr$nend[cbind(idx$g, idx$y)]
  X   <- do.call(rbind, lapply(xz, `[[`, "X")); colnames(X) <- dimnames(Theta)[[3]]
  list(X = X, z = unlist(lapply(xz, `[[`, "z")),
       id = data.frame(genotype = rep(tr$genotypes[idx$g], n), year = rep(tr$years[idx$y], n), day = sequence(n) - 1))
}

## A = X'X and b = X'z for every genotype and season. A linear fit needs nothing else,
## and cross-validation reduces to adding and subtracting these small matrices.
gram_stats <- function(Theta, tr) {
  G <- dim(tr$Y)[1]; ny <- length(tr$years); K <- dim(Theta)[3]
  A <- array(0, c(G, ny, K, K)); b <- array(0, c(G, ny, K))
  for (g in seq_len(G)) for (y in seq_len(ny)) {
    s <- trajectory_xz(Theta, tr, g, y)
    A[g, y, , ] <- crossprod(s$X); b[g, y, ] <- crossprod(s$X, s$z)
  }
  list(A = A, b = b, n = tr$nend, names = dimnames(Theta)[[3]])
}

## Sum the statistics over a set of seasons.
sum_years <- function(gram, years) {
  list(A = apply(gram$A[, years, , , drop = FALSE], c(1, 3, 4), sum),     # genotype x K x K
       b = apply(gram$b[, years, , drop = FALSE], c(1, 3), sum),           # genotype x K
       n = rowSums(gram$n[, years, drop = FALSE]))                         # genotype
}

## Scale the library columns to unit root-mean-square (coefficients become comparable).
scale_gram <- function(A, b, n, s = sqrt(diag(A) / n)) {
  list(A = A / n / outer(s, s), b = b / n / s, s = s)
}

## ---- Fitting methods: fit(A, b, ...) -> list(coef, support) ----------------------------
fit_ols <- function(A, b, ...) {
  list(coef = solve(A, b), support = seq_along(b))
}
fit_ridge <- function(A, b, lambda = 0.01, ...) {
  list(coef = solve(A + lambda * mean(diag(A)) * diag(nrow(A)), b), support = seq_along(b))
}
## Sequential thresholding: refit, drop the smallest coefficient, repeat until k terms remain.
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

## Lasso by coordinate descent on the scaled normal equations (soft thresholding).
fit_lasso <- function(A, b, lambda = 2e-4, iters = 2000, ...) {
  xi <- numeric(length(b))
  for (it in seq_len(iters)) {
    old <- xi
    for (j in seq_along(b)) {
      r <- b[j] - sum(A[j, -j] * xi[-j])
      xi[j] <- sign(r) * max(abs(r) - lambda, 0) / A[j, j]
    }
    if (max(abs(xi - old)) < 1e-9) break
  }
  list(coef = xi, support = which(xi != 0))
}

## ---- Genotypes: xi_g = mu + delta_g on a common support -------------------------------
## mu: ridge fit of the pooled data; delta_g: ridge-shrunk fit of the residual of genotype g.
fit_hier <- function(Ag, bg, support, lam_mu = 1e-4, lam_g = 0.01) {
  G <- dim(Ag)[1]; K <- dim(Ag)[2]; k <- length(support)
  Ap <- apply(Ag[, support, support, drop = FALSE], c(2, 3), sum); bp <- colSums(bg[, support, drop = FALSE])
  mu <- numeric(K); mu[support] <- solve(Ap + lam_mu * mean(diag(Ap)) * diag(k), bp)
  delta <- matrix(0, G, K)
  for (g in seq_len(G)) {
    Ags <- Ag[g, support, support]
    delta[g, support] <- solve(Ags + lam_g * mean(diag(Ags)) * diag(k), bg[g, support] - Ags %*% mu[support])
  }
  list(mu = mu, delta = delta, xi = sweep(delta, 2, mu, "+"), support = support)
}

## Fit on a set of seasons: `fitter` selects the shared support on the pooled problem,
## fit_hier gives mu and the deviations. Coefficients are returned in library units.
fit_years <- function(gram, years, fitter = fit_ssr, k = 12, lam_mu = 1e-4, lam_g = 0.01, ...) {
  st <- sum_years(gram, years); n <- sum(st$n)
  Apool <- apply(st$A, c(2, 3), sum); bpool <- colSums(st$b); s <- sqrt(diag(Apool) / n)
  Ag <- sweep(sweep(st$A, 2, s, "/"), 3, s, "/") / n
  bg <- sweep(st$b, 2, s, "/") / n
  support <- fitter(Apool / n / outer(s, s), bpool / n / s, k = k, ...)$support
  h <- fit_hier(Ag, bg, support, lam_mu, lam_g)
  list(xi = sweep(h$xi, 2, s, "/"), mu = h$mu / s, support = support, names = gram$names)
}

## ---- Simulation: for a weather-only library the solution is a cumulative sum -------------
## xi: genotype x term; Theta_y: day x term (one season); y0: initial biomass per genotype.
simulate <- function(xi, Theta_y, y0) {
  increments <- Theta_y[-nrow(Theta_y), , drop = FALSE] %*% t(xi)        # (days - 1) x genotype
  sweep(rbind(0, apply(increments, 2, cumsum)), 2, y0, "+")              # day x genotype
}

## ---- Validation: leave-one-year-out ---------------------------------------------------------
## Each season is predicted by a model fitted to the others, starting from their mean initial biomass.
cv_blocked <- function(Theta, tr, gram, fitter = fit_ssr, k = 12, lam_mu = 1e-4, lam_g = 0.01, ...) {
  ny <- length(tr$years); P <- array(NA_real_, c(dim(tr$Y)[1], ny, dim(Theta)[2])); support <- vector("list", ny)
  for (v in seq_len(ny)) {
    train <- setdiff(seq_len(ny), v)
    f <- fit_years(gram, train, fitter, k, lam_mu, lam_g, ...)
    P[, v, ] <- t(simulate(f$xi, Theta[v, , ], rowMeans(tr$Y[, train, 1])))
    support[[v]] <- f$support
  }
  list(pred = P, support = support)
}

## RMSE over whole trajectories, on a fixed day (array index; 99 is day 113 after sowing), and at maturity.
evaluate <- function(P, tr, fixed_day = 99) {
  se <- 0; cnt <- 0; final <- numeric(0)
  for (g in seq_len(dim(P)[1])) for (y in seq_len(dim(P)[2])) {
    n <- tr$nend[g, y]; d <- 2:(n + 1)                         # index 1 is the given initial value
    se <- se + sum((P[g, y, d] - tr$Y[g, y, d])^2); cnt <- cnt + n
    final <- c(final, P[g, y, n + 1] - tr$Y[g, y, n + 1])
  }
  c(trajectory_rmse = sqrt(se / cnt),
    fixed_day_rmse  = sqrt(mean((P[, , fixed_day] - tr$Y[, , fixed_day])^2, na.rm = TRUE)),
    maturity_rmse   = sqrt(mean(final^2)))
}

## RMSE on each day, over all genotype-season pairs that are still growing on that day.
rmse_by_day <- function(P, tr) {
  sapply(2:dim(P)[3], function(d) sqrt(mean((P[, , d] - tr$Y[, , d])^2, na.rm = TRUE)))
}

## Baseline without weather: each trajectory is the mean curve of the same genotype in the other
## seasons (held at its last value once all other seasons have ended), with the true initial value.
lookup_trajectory <- function(tr) {
  P <- array(NA_real_, dim(tr$Y))
  for (y in seq_along(tr$years)) {
    m <- apply(tr$Y[, -y, , drop = FALSE], c(1, 3), mean, na.rm = TRUE)
    P[, y, ] <- t(apply(m, 1, function(v) replace(v, is.nan(v), tail(v[!is.nan(v)], 1))))
    P[, y, 1] <- tr$Y[, y, 1]
  }
  P
}

## ---- Interpretation: accumulated local effects (ALE) --------------------------------------
## Standardised drivers on the days that occur in the data (rows = season x day).
observed_points <- function(drv, tr) {
  do.call(rbind, lapply(seq_along(drv$years), function(y) {
    n <- max(tr$nend[, y]); matrix(drv$U[y, 1:n, ], n, dimnames = list(NULL, drv$drivers))
  }))
}
## Effect of driver j on the growth rate. Within each of K quantile intervals, driver j is moved
## from the lower to the upper end while the other drivers keep their observed values; the mean
## changes in the rate are accumulated and centred. rate: function(Z) -> growth rate.
ale_curve <- function(rate, Z, j, K = 24) {
  edge <- unique(quantile(Z[, j], seq(0, 1, length.out = K + 1))); nb <- length(edge) - 1
  bin  <- pmin(pmax(findInterval(Z[, j], edge), 1), nb)
  Zu <- Zl <- Z; Zu[, j] <- edge[bin + 1]; Zl[, j] <- edge[bin]
  count  <- tabulate(bin, nb)
  change <- vapply(split(rate(Zu) - rate(Zl), factor(bin, levels = seq_len(nb))), sum, 0) / pmax(count, 1)
  A <- c(0, cumsum(change)); mid <- (A[-1] + A[-length(A)]) / 2
  data.frame(z = edge, effect = A - sum(mid * count) / sum(count))
}

## ---- Genotype ranking ---------------------------------------------------------------------
## Pd, Yd: genotype x season matrices of predicted and APSIM biomass on one day. For each season:
## Spearman correlation of the two rankings and overlap of the predicted and true top 10 %.
rank_scores <- function(Pd, Yd, top = 0.1) {
  k <- max(1, round(top * nrow(Yd)))
  do.call(rbind, lapply(seq_len(ncol(Yd)), function(y) data.frame(
    year = y, rho = cor(Pd[, y], Yd[, y], method = "spearman"),
    top = length(intersect(order(-Pd[, y])[1:k], order(-Yd[, y])[1:k])) / k)))
}
## Baseline without weather: the genotype's mean on that day in the other seasons.
lookup_baseline <- function(Yd) (rowSums(Yd) - Yd) / (ncol(Yd) - 1)
