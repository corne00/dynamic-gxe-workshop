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

## path(X, dy, values) -> list of coefficient vectors, one per value
path_stlsq <- function(X, dy, values) {
  Z <- cbind(`1` = 1, X); G <- crossprod(Z); rhs <- drop(crossprod(Z, dy))
  size <- col_sd(X) / sd(dy)
  lapply(values, function(threshold) {
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
  })
}

path_ridge <- function(X, dy, values) {
  n <- nrow(X); mx <- colMeans(X)
  var_x <- col_sd(X)^2 * (n - 1) / n
  sy <- sd(dy) * sqrt((n - 1) / n)
  s <- 1 / sqrt(var_x)
  e <- eigen(s * (crossprod(X) / n - tcrossprod(mx)) * rep(s, each = length(s)),
             symmetric = TRUE)
  r <- drop(s * (crossprod(X, dy) / n - mx * mean(dy)))
  q <- crossprod(e$vectors, r)
  lapply(values, function(lambda) {
    b <- s * drop(e$vectors %*% (q / (e$values + lambda / sy)))
    c(`1` = mean(dy) - sum(mx * b), setNames(b, colnames(X)))
  })
}

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
