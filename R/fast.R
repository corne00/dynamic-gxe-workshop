## Faster routes to the numbers behind the precomputed results. They return the
## same values as cross_validate() with the fit_* functions of R/functions.R:
## genotypes are spread over cores, and a path over a hyperparameter is fitted
## in one call (one glmnet path, one eigen decomposition, one leaps search).

cross_validate_par <- function(Theta, tr, fit, genotypes = tr$genotypes,
                               cores = parallel::detectCores()) {
  chunks <- split(genotypes, rep_len(seq_len(cores), length(genotypes)))
  parts <- parallel::mclapply(chunks, function(g) cross_validate(Theta, tr, fit, g),
                              mc.cores = cores)
  P <- array(NA_real_, dim(tr$Y))
  for (p in parts) { i <- !is.na(p[, 1, 1]); P[i, , ] <- p[i, , ] }
  P
}

## Moment-based fits. G = crossprod(cbind(1, X, dy)) holds everything that
## OLS, ridge and STLSQ need, and the moments of a training set are the totals
## minus those of the held-out year, so no fold touches the data again.
## An mpath(G, values) returns one coefficient vector (intercept first) per value.
moments <- function(X, dy) {
  G <- crossprod(cbind(1, X, dy))
  dimnames(G) <- list(NULL, NULL)
  structure(G, names_x = colnames(X))
}
moment_parts <- function(G) {
  K <- ncol(G) - 2; n <- G[1, 1]
  list(K = K, n = n, zz = G[1:(K + 1), 1:(K + 1), drop = FALSE],
       zy = G[1:(K + 1), K + 2], sum_x = G[1, 2:(K + 1)], sum_y = G[1, K + 2],
       ss_y = G[K + 2, K + 2],
       var_x = (diag(G)[2:(K + 1)] - G[1, 2:(K + 1)]^2 / n) / (n - 1))
}
labelled <- function(b, names_x) setNames(b, c("1", names_x))

mpath_ols <- function(G, values = 0) {
  m <- moment_parts(G); names_x <- attr(G, "names_x")
  rep(list(labelled(drop(solve(m$zz, m$zy)), names_x)), length(values))
}

mpath_stlsq <- function(G, values) {
  m <- moment_parts(G); names_x <- attr(G, "names_x")
  sd_y <- sqrt((m$ss_y - m$sum_y^2 / m$n) / (m$n - 1))
  size <- sqrt(m$var_x) / sd_y
  idx <- seq_len(m$K + 1)
  lapply(values, function(threshold) {
    active <- seq_len(m$K)                       # columns of X
    repeat {
      i <- c(1, active + 1)
      b <- drop(solve(m$zz[i, i], m$zy[i]))
      keep <- abs(b[-1]) * size[active] >= threshold
      if (all(keep) || !any(keep)) break
      active <- active[keep]
    }
    full <- numeric(m$K + 1); full[i] <- b
    labelled(full, names_x)
  })
}

mpath_ridge <- function(G, values) {
  m <- moment_parts(G); names_x <- attr(G, "names_x"); n <- m$n
  mx <- m$sum_x / n; my <- m$sum_y / n
  var_x <- m$var_x * (n - 1) / n
  sy <- sqrt((m$ss_y - m$sum_y^2 / n) / n)       # population sd of dy
  s <- 1 / sqrt(var_x)
  C <- m$zz[-1, -1, drop = FALSE] / n - tcrossprod(mx)
  e <- eigen(s * C * rep(s, each = m$K), symmetric = TRUE)
  r <- drop(s * (m$zy[-1] / n - mx * my))
  q <- crossprod(e$vectors, r)
  lapply(values, function(lambda) {
    b <- s * drop(e$vectors %*% (q / (e$values + lambda / sy)))
    labelled(c(my - sum(mx * b), b), names_x)
  })
}

## The X form of the same paths.
path_stlsq <- function(X, dy, values) mpath_stlsq(moments(X, dy), values)
path_ridge <- function(X, dy, values) mpath_ridge(moments(X, dy), values)

path_glmnet <- function(alpha) function(X, dy, values) {
  o <- order(values, decreasing = TRUE)
  fit <- glmnet::glmnet(X, dy, alpha = alpha, lambda = values[o],
                        type.gaussian = "covariance")
  B <- as.matrix(coef(fit))
  stopifnot(ncol(B) == length(values))
  lapply(match(seq_along(values), o), function(i) named(B[, i], X))
}

path_subset <- function(method) function(X, dy, values) {
  reg <- leaps::regsubsets(X, dy, nvmax = max(values), method = method)
  lapply(values, function(k) {
    b <- coef(reg, k)
    full <- setNames(numeric(ncol(X) + 1), c("1", colnames(X)))
    full[c("1", names(b)[-1])] <- b
    full
  })
}

## Trajectory RMSE of leave-one-year-out validation for every value of a path;
## equals evaluate(cross_validate(...))[["trajectory_rmse"]] per value.
cv_rmse_path <- function(Theta, tr, path, values, genotypes,
                         cores = parallel::detectCores()) {
  S <- length(tr$years)
  one <- function(gs) {
    ss <- numeric(length(values)); cnt <- 0
    for (g in match(gs, tr$genotypes)) {
      d <- make_design(Theta, tr, tr$genotypes[g])
      folds <- lapply(seq_len(S), function(v) {
        train <- d$year != v
        path(d$X[train, ], d$dy[train], values)
      })
      y0 <- (sum(tr$Y[g, , 1]) - tr$Y[g, , 1]) / (S - 1)
      for (i in seq_along(values)) {
        P <- t(simulate(t(sapply(folds, `[[`, i)), Theta, y0))
        e <- (P - tr$Y[g, , ])[, -1]
        ss[i] <- ss[i] + sum(e^2, na.rm = TRUE)
        if (i == 1) cnt <- cnt + sum(!is.na(e))
      }
    }
    c(ss, cnt)
  }
  chunks <- split(genotypes, rep_len(seq_len(cores), length(genotypes)))
  r <- Reduce(`+`, parallel::mclapply(chunks, one, mc.cores = cores))
  sqrt(r[seq_along(values)] / r[length(values) + 1])
}

## Leave-one-year-out validation with a moment-based path. Returns the
## trajectory RMSE per value (what = "rmse"), or the held-out predictions
## P[genotype, year, day] per value (what = "P"); equals the X-form results.
cv_moments <- function(Theta, tr, mpath, values, genotypes = tr$genotypes,
                       what = "rmse", cores = parallel::detectCores()) {
  S <- length(tr$years); nv <- length(values)
  one <- function(gs) {
    ss <- numeric(nv); cnt <- 0
    P <- if (what == "P") replicate(nv, array(NA_real_, dim(tr$Y)), FALSE)
    for (g in match(gs, tr$genotypes)) {
      d <- make_design(Theta, tr, tr$genotypes[g])
      by_year <- lapply(seq_len(S), function(v) {
        i <- d$year == v
        crossprod(cbind(1, d$X[i, , drop = FALSE], d$dy[i]))
      })
      total <- Reduce(`+`, by_year)
      folds <- lapply(by_year, function(Gv) {
        G <- total - Gv; attr(G, "names_x") <- colnames(d$X)
        mpath(G, values)
      })
      y0 <- (sum(tr$Y[g, , 1]) - tr$Y[g, , 1]) / (S - 1)
      for (i in seq_len(nv)) {
        Pi <- t(simulate(t(sapply(folds, `[[`, i)), Theta, y0))
        if (what == "P") { P[[i]][g, , ] <- Pi; next }
        e <- (Pi - tr$Y[g, , ])[, -1]
        ss[i] <- ss[i] + sum(e^2, na.rm = TRUE)
        if (i == 1) cnt <- cnt + sum(!is.na(e))
      }
    }
    if (what == "P") P else c(ss, cnt)
  }
  chunks <- split(genotypes, rep_len(seq_len(cores), length(genotypes)))
  parts <- parallel::mclapply(chunks, one, mc.cores = cores)
  if (what == "rmse") {
    r <- Reduce(`+`, parts)
    return(sqrt(r[seq_len(nv)] / r[nv + 1]))
  }
  lapply(seq_len(nv), function(i) {
    P <- array(NA_real_, dim(tr$Y))
    for (p in parts) { k <- !is.na(p[[i]][, 1, 1]); P[k, , ] <- p[[i]][k, , ] }
    P
  })
}
