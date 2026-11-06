## Precomputed results for the interactive figures (written to
## data/interactive/*.json). Run from the repository root; takes a few minutes.
library(jsonlite)
source("R/pipeline.R")

out <- function(x, name, digits = 4) {
  write_json(x, file.path("data/interactive", name), digits = digits,
             auto_unbox = TRUE, na = "null")
}
rounded <- function(a, digits = 2) round(a, digits)
three <- spread_genotypes(tr)
some <- union(three, tr$genotypes[seq(1, 199, by = 5)])

## 1. Hyperparameters: scaled coefficients per genotype and held-out error ----
## The error curve is the mean over cv_genotypes; the coefficients are shown
## for one genotype at a time, and the terms kept differ between genotypes.
cv_genotypes <- tr$genotypes[seq(1, 199, by = 10)]
designs <- lapply(cv_genotypes, function(g) make_design(Theta, tr, g))
setting <- function(fit) {
  fits <- lapply(designs, function(d) {
    b <- fit(d$X, d$dy)
    list(effect = unname(b[-1] * apply(d$X, 2, sd) / sd(d$dy)),
         nterms = sum(b[-1] != 0))
  })
  list(effect = lapply(fits, `[[`, "effect"),
       nterms = sapply(fits, `[[`, "nterms"),
       rmse = evaluate(cross_validate(Theta, tr, fit, cv_genotypes),
                       tr)[["trajectory_rmse"]])
}
grid <- list(
  STLSQ = list(param = "threshold", values = c(0.005, 0.01, 0.02, 0.03, 0.05,
                                                0.08, 0.12, 0.2, 0.3),
               fit = function(v) function(X, dy) fit_stlsq(X, dy, v)),
  lasso = list(param = "lambda", values = 10^seq(-5, -2.5, length.out = 11),
               fit = function(v) function(X, dy) fit_lasso(X, dy, v)),
  `elastic net` = list(param = "lambda",
                       values = 10^seq(-5, -2.5, length.out = 11),
                       fit = function(v) function(X, dy) fit_enet(X, dy, v)),
  ridge = list(param = "lambda", values = 10^seq(-5, -1, length.out = 11),
               fit = function(v) function(X, dy) fit_ridge(X, dy, v)),
  `best subset` = list(param = "terms", values = 1:12,
                       fit = function(v) function(X, dy) fit_subset(X, dy, v)),
  backward = list(param = "terms", values = 1:12,
                  fit = function(v) function(X, dy) {
                    fit_subset(X, dy, v, "backward")
                  }))
methods <- lapply(grid, function(m) {
  s <- lapply(m$values, function(v) setting(m$fit(v)))
  list(param = m$param, values = m$values,
       effect = lapply(s, `[[`, "effect"),
       nterms = t(sapply(s, `[[`, "nterms")), rmse = sapply(s, `[[`, "rmse"))
})
out(list(terms = colnames(designs[[1]]$X), genotypes = cv_genotypes,
         methods = methods), "hyper.json", 5)

## 2. Held-out trajectories of the sparse model and of the baseline -----------
P <- readRDS("precomputed/cv_stlsq.rds")
base <- lookup_trajectory(tr)
pick <- match(some, tr$genotypes)
pack <- function(A) {
  lapply(seq_along(pick), function(i) lapply(seq_along(tr$years), function(s) {
    n <- tr$nend[pick[i], s] + 1
    rounded(A[pick[i], s, 1:n])
  }))
}
out(list(genotypes = some, years = tr$years, das0 = tr$das0,
         three = three, observed = pack(tr$Y), model = pack(P),
         baseline = pack(base)), "trajectories.json")

## 3. ALE curves of all genotypes ---------------------------------------------
ale <- readRDS("precomputed/ale_genotypes.rds")
ale_list <- lapply(split(ale, ale$driver), function(a) {
  lapply(split(a, a$genotype), function(z) {
    list(x = signif(z$x, 5), effect = signif(z$effect, 4))
  })
})
out(list(drivers = names(ale_list), curves = ale_list, three = three),
    "ale.json")

## 4. Library and coefficients (the page computes the contributions) ----------
coefs <- sapply(tr$genotypes, function(g) {
  dd <- make_design(Theta, tr, g)
  fit_stlsq(dd$X, dd$dy)
})
colnames(coefs) <- tr$genotypes
out(list(terms = dimnames(Theta)[[3]], drivers = drv$drivers, years = tr$years,
         das0 = tr$das0, genotypes = tr$genotypes, three = three,
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
         observed = pack_day(tr$Y), model = pack_day(P),
         baseline = pack_day(base)), "ranking.json")
