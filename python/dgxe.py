"""Sparse identification of crop growth: Python counterpart of R/functions.R.

  weather --make_drivers--> drv --make_library--> Theta   (candidate terms)
  biomass --make_trajectories--> tr
  Theta, tr --make_design--> X, dy                (regression of increments)
  fit_*(X, dy) --> coefficients --simulate--> prediction --evaluate--> RMSE

numpy, pandas, scipy and scikit-learn. Indices are 0-based; days are given
as days after sowing (das).
"""
from itertools import combinations
from types import SimpleNamespace

import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet

DRIVERS = ["TT", "DL", "RAW", "RAD"]
GREY, BLUE, CYAN, BURGUNDY = "#D1D5DB", "#014E81", "#00A6D6", "#8B1E3F"


# ---- Drivers: array U[season, day, driver], standardised -------------------
def make_drivers(wx, drivers=DRIVERS):
    years = np.sort(wx["year"].unique())
    das0 = wx["das"].min()
    U = np.full((len(years), wx["das"].max() - das0 + 1, len(drivers)), np.nan)
    y = np.searchsorted(years, wx["year"])
    U[y, (wx["das"] - das0).to_numpy(), :] = wx[drivers].to_numpy()
    for s in range(len(years)):                  # seasons end on different days
        for j in range(len(drivers)):            # pad with the last value
            u = U[s, :, j]
            u[np.isnan(u)] = u[~np.isnan(u)][-1]
    center = U.mean(axis=(0, 1))
    scale = U.std(axis=(0, 1))
    return SimpleNamespace(U=(U - center) / scale, years=years, das0=das0,
                           drivers=list(drivers), center=center, scale=scale)


# ---- Library of candidate terms --------------------------------------------
# Z: n x m matrix of standardised drivers -> n x K matrix of terms. Powers
# 1..degree of every driver; pairs = True adds products of pairs of drivers.
# The constant term is the intercept of the regression and is not included.
def library_matrix(Z, names, degree=3, pairs=False):
    def term(x, p):
        return x if p == 1 else f"{x}^{p}"
    cols = {}
    for j, d in enumerate(names):
        for p in range(1, degree + 1):
            cols[term(d, p)] = Z[:, j] ** p
    if pairs:
        for i, j in combinations(range(len(names)), 2):
            for p in (1, 2):
                for q in (1, 2):
                    name = f"{term(names[i], p)}*{term(names[j], q)}"
                    cols[name] = Z[:, i] ** p * Z[:, j] ** q
    return np.column_stack(list(cols.values())), list(cols)


# Theta[season, day, term]
def make_library(drv, degree=3, pairs=False):
    n_season, n_day, n_driver = drv.U.shape
    Z = drv.U.reshape(-1, n_driver)
    L, terms = library_matrix(Z, drv.drivers, degree, pairs)
    return SimpleNamespace(array=L.reshape(n_season, n_day, -1), terms=terms)


# ---- Trajectories: array Y[genotype, season, day], NaN after maturity ------
# nend[g, s] is the number of daily steps of trajectory (g, s).
def make_trajectories(traj, drv, var="biomass"):
    genotypes = np.sort(traj["genotype"].unique())
    Y = np.full((len(genotypes), len(drv.years), drv.U.shape[1]), np.nan)
    g = np.searchsorted(genotypes, traj["genotype"])
    s = np.searchsorted(drv.years, traj["year"])
    Y[g, s, traj["das"] - drv.das0] = traj[var]
    return SimpleNamespace(Y=Y, nend=(~np.isnan(Y)).sum(axis=2) - 1,
                           genotypes=genotypes, years=drv.years,
                           das0=drv.das0)


# ---- Regression problem for one genotype -----------------------------------
# Rows: all days of all seasons. X: library at day t. dy: y(t + 1) - y(t).
def make_design(Theta, tr, genotype):
    g = list(tr.genotypes).index(genotype)
    n = tr.nend[g]
    seasons = range(len(tr.years))
    X = np.vstack([Theta.array[s, :n[s], :] for s in seasons])
    dy = np.concatenate([np.diff(tr.Y[g, s, :n[s] + 1]) for s in seasons])
    season = np.repeat(np.arange(len(tr.years)), n)
    return SimpleNamespace(X=X, dy=dy, season=season, terms=Theta.terms)


# ---- Fitting methods: fit(X, dy) -> coefficients (intercept first) ---------
def fit_ols(X, dy):
    return np.linalg.lstsq(np.column_stack([np.ones(len(dy)), X]), dy,
                           rcond=None)[0]


# Ridge (alpha = 0), lasso (alpha = 1) and elastic net with the penalty of
# glmnet, lam * ((1 - alpha) / 2 * b'b / sd(dy) + alpha * |b|_1), on terms
# scaled to unit variance (glmnet scales the squared penalty by sd(dy)).
def fit_glmnet(X, dy, alpha, lam):
    sx = X.std(axis=0)
    l1, l2 = lam * alpha, lam * (1 - alpha) / dy.std()
    model = ElasticNet(alpha=l1 + l2, l1_ratio=l1 / (l1 + l2), tol=1e-12,
                       max_iter=100000)
    model.fit((X - X.mean(axis=0)) / sx, dy - dy.mean())
    b = model.coef_ / sx
    return np.concatenate([[dy.mean() - X.mean(axis=0) @ b], b])


def fit_ridge(X, dy, lam=1e-3): return fit_glmnet(X, dy, 0, lam)
def fit_lasso(X, dy, lam=3e-4): return fit_glmnet(X, dy, 1, lam)
def fit_enet(X, dy, lam=6e-4): return fit_glmnet(X, dy, 0.5, lam)


# Best subset (exhaustive) or backward elimination with k terms. Residual sums
# of squares come from the centred Gram matrix, so no refit is needed.
def fit_subset(X, dy, k=8, method="exhaustive"):
    Xc, yc = X - X.mean(axis=0), dy - dy.mean()
    G, c = Xc.T @ Xc, Xc.T @ yc

    def rss(s):
        s = list(s)
        return yc @ yc - c[s] @ np.linalg.solve(G[np.ix_(s, s)], c[s])

    if method == "exhaustive":
        keep = min(combinations(range(X.shape[1]), k), key=rss)
    else:
        keep = list(range(X.shape[1]))
        while len(keep) > k:
            keep.remove(min(keep, key=lambda j: rss([i for i in keep
                                                     if i != j])))
    keep = list(keep)
    b = np.zeros(X.shape[1])
    b[keep] = np.linalg.solve(G[np.ix_(keep, keep)], c[keep])
    return np.concatenate([[dy.mean() - X.mean(axis=0) @ b], b])


# Sequential thresholded least squares (the original SINDy algorithm): refit
# on the active terms and remove those whose scaled coefficient is below the
# threshold, until no term is removed.
def fit_stlsq(X, dy, threshold=0.05):
    size = X.std(axis=0, ddof=1) / dy.std(ddof=1)    # common scale
    active = np.arange(X.shape[1])
    while True:
        b = fit_ols(X[:, active], dy)
        keep = np.abs(b[1:]) * size[active] >= threshold
        if keep.all() or not keep.any():
            break
        active = active[keep]
    full = np.zeros(X.shape[1] + 1)
    full[0], full[1 + active] = b[0], b[1:]
    return full


# Methods compared in the workshop (hyperparameters fixed for simplicity).
fit_methods = {
    "OLS": fit_ols,
    "ridge": fit_ridge,
    "lasso": fit_lasso,
    "elastic net": fit_enet,
    "best subset (k = 8)": lambda X, dy: fit_subset(X, dy, 8),
    "backward (k = 8)": lambda X, dy: fit_subset(X, dy, 8, "backward"),
    "STLSQ": fit_stlsq,
}


# ---- Simulation: for a weather-only library the solution is a cumulative sum
# coef: intercept and term coefficients; Theta_s: day x term (one season).
def simulate(coef, Theta_s, y0):
    rate = coef[0] + Theta_s[:-1] @ coef[1:]
    return y0 + np.concatenate([[0], np.cumsum(rate)])


# Genotypes with the lowest, median and highest mean final biomass.
def spread_genotypes(tr):
    g, s = np.indices(tr.nend.shape)
    final = tr.Y[g, s, tr.nend].mean(axis=1)
    ranked = tr.genotypes[np.argsort(final)]
    return ranked[[0, int(np.ceil(len(ranked) / 2)) - 1, -1]]


# ---- Validation: leave-one-season-out --------------------------------------
# Each season is predicted by a model fitted to the other seasons of the same
# genotype, starting from the mean initial biomass of those seasons.
# Returns P[genotype, season, day].
def cross_validate(Theta, tr, fit, genotypes=None):
    P = np.full(tr.Y.shape, np.nan)
    for g in (tr.genotypes if genotypes is None else genotypes):
        i = list(tr.genotypes).index(g)
        d = make_design(Theta, tr, g)
        for s in range(len(tr.years)):
            train = d.season != s
            coef = fit(d.X[train], d.dy[train])
            y0 = np.delete(tr.Y[i, :, 0], s).mean()
            P[i, s] = simulate(coef, Theta.array[s], y0)
    return P


# RMSE over whole trajectories, on a fixed day (days after sowing) and at
# maturity. Genotypes that were not predicted are ignored.
def evaluate(P, tr, fixed_day=113):
    err = P - tr.Y
    err[:, :, 0] = np.nan                        # the initial value is given
    g, s = np.indices(tr.nend.shape)

    def rmse(e): return np.sqrt(np.nanmean(e ** 2))
    return {"trajectory_rmse": rmse(err),
            "fixed_day_rmse": rmse(err[:, :, fixed_day - tr.das0]),
            "maturity_rmse": rmse(err[g, s, tr.nend])}


# RMSE on each day, over all genotype-season pairs still growing on that day.
def rmse_by_day(P, tr):
    return np.sqrt(np.nanmean((P - tr.Y) ** 2, axis=(0, 1)))[1:]


# Baseline without weather: each trajectory is the mean curve of the same
# genotype in the other seasons (held at its last value once all other
# seasons have ended), with the true initial value.
def lookup_trajectory(tr):
    P = np.full(tr.Y.shape, np.nan)
    for s in range(len(tr.years)):
        others = np.delete(tr.Y, s, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = np.nansum(others, axis=1) / (~np.isnan(others)).sum(axis=1)
        m[(~np.isnan(others)).sum(axis=1) == 0] = np.nan
        for g in range(m.shape[0]):
            ok = ~np.isnan(m[g])
            m[g, ~ok] = m[g, ok][-1]
        P[:, s] = m
        P[:, s, 0] = tr.Y[:, s, 0]
    return P


# ---- Genotype ranking -------------------------------------------------------
# Pd, Yd: genotype x season matrices of predicted and APSIM biomass on one
# day. For each season: Spearman correlation of the two rankings and overlap
# of the predicted and true top 10 %.
def rank_scores(Pd, Yd, top=0.1):
    from scipy.stats import spearmanr
    k = max(1, round(top * Yd.shape[0]))
    rho, overlap = [], []
    for s in range(Yd.shape[1]):
        rho.append(spearmanr(Pd[:, s], Yd[:, s])[0])
        best = set(np.argsort(-Pd[:, s])[:k]) & set(np.argsort(-Yd[:, s])[:k])
        overlap.append(len(best) / k)
    return np.array(rho), np.array(overlap)


# Baseline without weather: the genotype's mean on that day in the other
# seasons.
def lookup_baseline(Yd):
    return (Yd.sum(axis=1, keepdims=True) - Yd) / (Yd.shape[1] - 1)


# ---- Interpretation --------------------------------------------------------
# Standardised drivers on the days that occur in the data (season x day rows).
def observed_points(drv, tr):
    last = tr.nend.max(axis=0)
    rows = [drv.U[s, :last[s]] for s in range(len(tr.years))]
    return pd.DataFrame(np.vstack(rows), columns=drv.drivers)


# Growth rate of one genotype as a function of the drivers (a data frame).
class GrowthRate:
    def __init__(self, coef, pairs=False):
        self.coef, self.pairs = coef, pairs

    def predict(self, Z):
        L, _ = library_matrix(Z.to_numpy(), list(Z.columns), pairs=self.pairs)
        return self.coef[0] + L @ self.coef[1:]
