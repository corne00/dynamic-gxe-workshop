## Precomputed results for the interactive figures (written to
## data/interactive/*.json). Run from the repository root; takes a few minutes.
library(jsonlite)
source("R/pipeline.R")
source("scripts/helpers.R")
source("R/fast.R")

out <- function(x, name, digits = 4) {
  write_json(x, file.path("data/interactive", name), digits = digits,
             auto_unbox = TRUE, na = "null")
}
rounded <- function(a, digits = 2) round(a, digits)
selected <- spread_genotypes(tr)
some <- union(selected, tr$genotypes[seq(1, 199, by = 5)])

## 1. Hyperparameters: scaled coefficients per genotype and held-out error ----
## The error curve is the mean over cv_genotypes; the coefficients are shown
## for one genotype at a time, and the terms kept differ between genotypes.
cv_genotypes <- tr$genotypes[seq(1, 199, by = 10)]
designs <- lapply(cv_genotypes, function(g) make_design(Theta, tr, g))
sd_x <- lapply(designs, function(d) apply(d$Theta, 2, sd))
sd_y <- sapply(designs, function(d) sd(d$dy))
grid <- list(
  STLSQ = list(param = "threshold", path = path_stlsq, mpath = mpath_stlsq,
               values = signif(10^seq(log10(0.004), log10(0.5), length.out = 30), 2)),
  lasso = list(param = "lambda", path = path_glmnet(1),
               values = 10^seq(-5, -0.9, length.out = 30)),
  `elastic net` = list(param = "lambda", path = path_glmnet(0.5),
                       values = 10^seq(-5, -0.9, length.out = 30)),
  ridge = list(param = "lambda", path = path_ridge, mpath = mpath_ridge,
               values = 10^seq(-5, 0, length.out = 30)),
  `best subset` = list(param = "terms", path = path_subset("exhaustive"),
                       values = 1:24),
  backward = list(param = "terms", path = path_subset("backward"),
                  values = 1:24))
## Held-out trajectories (leave-one-year-out) of a few genotypes and years for
## every value of a method, every second day, to show next to the error curve.
contrast <- c("g181", "g171", "g021", "g141", "g151")
traj_years <- c(1984, 1998, 2008)
step <- 2
year_index <- match(traj_years, tr$years)
days_kept <- function(g, s) seq(1, tr$nend[g, s] + 1, by = step)
held_out <- function(m, name) {
  g <- match(name, tr$genotypes); S <- length(tr$years)
  d <- make_design(Theta, tr, name)
  folds <- lapply(seq_len(S), function(v) {
    train <- d$year != v
    m$path(d$Theta[train, ], d$dy[train], m$values)
  })
  y0 <- (sum(tr$Y[g, , 1]) - tr$Y[g, , 1]) / (S - 1)
  lapply(seq_along(m$values), function(i) {
    P <- t(simulate(t(sapply(folds, `[[`, i)), Theta, y0))
    lapply(year_index, function(s) {
      p <- P[s, days_kept(g, s)]
      p[!is.finite(p) | abs(p) > 50] <- NA
      round(p, 2)
    })
  })
}
methods <- lapply(grid, function(m) {
  fits <- lapply(designs, function(d) m$path(d$Theta, d$dy, m$values))
  held <- lapply(contrast, function(name) held_out(m, name))
  pred <- lapply(seq_along(m$values), function(i) lapply(held, `[[`, i))
  i_val <- seq_along(m$values); i_gen <- seq_along(designs)
  list(param = m$param, values = m$values, pred = pred,
       effect = lapply(i_val, function(i) lapply(i_gen, function(g) {
         unname(fits[[g]][[i]][-1] * sd_x[[g]] / sd_y[g])
       })),
       coef = lapply(i_val, function(i) lapply(i_gen, function(g) {
         signif(unname(fits[[g]][[i]][-1]), 4)
       })),
       intercept = lapply(i_val, function(i) sapply(i_gen, function(g) {
         signif(unname(fits[[g]][[i]][1]), 4)
       })),
       nterms = t(sapply(i_val, function(i) sapply(i_gen, function(g) {
         sum(fits[[g]][[i]][-1] != 0)
       }))),
       rmse = if (is.null(m$mpath)) cv_rmse_path(Theta, tr, m$path, m$values, cv_genotypes)
               else cv_moments(Theta, tr, m$mpath, m$values, cv_genotypes))
})
observed <- lapply(match(contrast, tr$genotypes), function(g) {
  lapply(year_index, function(s) round(tr$Y[g, s, days_kept(g, s)], 2))
})
out(list(terms = colnames(designs[[1]]$Theta), genotypes = cv_genotypes,
         methods = methods,
         trajectories = list(genotypes = contrast, years = traj_years,
                             das0 = tr$das0, step = step,
                             observed = observed)), "hyper.json", 5)

## 2. Held-out trajectories of the sparse model -------------------------------
P <- readRDS("precomputed/cv_stlsq.rds")
contrast <- spread_genotypes(tr, 5)
pick <- match(contrast, tr$genotypes)
pack <- function(A) {
  lapply(seq_along(pick), function(i) lapply(seq_along(tr$years), function(s) {
    n <- tr$nend[pick[i], s] + 1
    rounded(A[pick[i], s, 1:n])
  }))
}
out(list(genotypes = contrast, years = tr$years, das0 = tr$das0,
         observed = pack(tr$Y), model = pack(P)), "trajectories.json")

## 3. ALE curves of all genotypes ---------------------------------------------
ale <- readRDS("precomputed/ale_genotypes.rds")
ale_list <- lapply(split(ale, ale$driver), function(a) {
  lapply(split(a, a$genotype), function(z) {
    list(x = signif(z$x, 5), effect = signif(z$effect, 4))
  })
})
out(list(drivers = names(ale_list), curves = ale_list, selected = selected),
    "ale.json")

## 4. Library and coefficients (the page computes the contributions) ----------
coefs <- sapply(tr$genotypes, function(g) {
  dd <- make_design(Theta, tr, g)
  fit_stlsq(dd$Theta, dd$dy)
})
colnames(coefs) <- tr$genotypes
out(list(terms = dimnames(Theta)[[3]], drivers = drv$drivers, years = tr$years,
         das0 = tr$das0, genotypes = tr$genotypes, selected = selected,
         nend = lapply(seq_along(tr$genotypes), function(g) tr$nend[g, ]),
         theta = lapply(seq_along(tr$years), function(s) {
           lapply(1:dim(Theta)[2], function(i) round(Theta[s, i, ], 4))
         }),
         coefs = lapply(tr$genotypes, function(g) signif(coefs[, g], 5))),
    "contrib.json")

## 5. Ranking of the genotypes on a given day ----------------------------------
days <- seq(29, 127, by = 7)                    # days after sowing
index <- days - tr$das0 + 1
pack_day <- function(A) {
  lapply(seq_along(tr$years), function(s) lapply(index, function(i) {
    rounded(A[, s, i])
  }))
}
out(list(days = days, years = tr$years, genotypes = tr$genotypes,
         observed = pack_day(tr$Y), model = pack_day(P)), "ranking.json")
