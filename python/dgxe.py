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
from scipy.integrate import solve_ivp
from sklearn.linear_model import ElasticNet

DRIVERS = ["Tmin", "Tmax", "Rain", "RAD", "VP", "DL", "TT", "RAW"]
GREY, BLUE, CYAN, BURGUNDY = "#D1D5DB", "#014E81", "#00A6D6", "#8B1E3F"


# ---- Drivers: array U[year, day, driver], standardised -------------------
def make_drivers(wx, drivers=DRIVERS):
    years = np.sort(wx["year"].unique())
    das0 = wx["das"].min()
    U = np.full((len(years), wx["das"].max() - das0 + 1, len(drivers)), np.nan)
    y = np.searchsorted(years, wx["year"])
    U[y, (wx["das"] - das0).to_numpy(), :] = wx[drivers].to_numpy()
    for s in range(len(years)):                  # years end on different days
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


# Theta[year, day, term]
def make_library(drv, degree=3, pairs=False):
    n_year, n_day, n_driver = drv.U.shape
    Z = drv.U.reshape(-1, n_driver)
    L, terms = library_matrix(Z, drv.drivers, degree, pairs)
    return SimpleNamespace(array=L.reshape(n_year, n_day, -1), terms=terms)


# ---- Extra terms: functions of the raw weather, standardised like the drivers
# exprs: dict name -> expression in Python syntax (^ is accepted for powers).
# Available are the columns of wx, Tmean, VPD (hPa) and das, and the functions
# roll(x, k) (mean of the last k days), cum(x) and lag(x, k); max and min act
# elementwise.
def extra_terms(wx, exprs, years=None):
    years = sorted(wx["year"].unique()) if years is None else years
    das0 = wx["das"].min()
    nd = wx["das"].max() - das0 + 1

    def roll(x, k):
        run = np.concatenate([[0.0], np.cumsum(x)])
        t = np.arange(len(x))
        start = np.maximum(0, t - k + 1)
        return (run[t + 1] - run[start]) / (t + 1 - start)

    funs = dict(exp=np.exp, log=np.log, sqrt=np.sqrt, abs=np.abs, sin=np.sin,
                cos=np.cos, tan=np.tan, tanh=np.tanh, pi=np.pi, max=np.maximum,
                min=np.minimum, cum=np.cumsum, roll=roll,
                lag=lambda x, k: x[np.maximum(np.arange(len(x)) - k, 0)])
    X = np.full((len(years), nd, len(exprs)), np.nan)
    for s, year in enumerate(years):
        w = wx[wx["year"] == year]
        v = {c: np.concatenate([w[c], np.full(nd - len(w), w[c].iloc[-1])])
             for c in w.columns if c not in ("year", "das")}
        v["das"] = das0 + np.arange(nd)
        v["Tmean"] = (v["Tmin"] + v["Tmax"]) / 2
        v["VPD"] = 6.108 * np.exp(17.27 * v["Tmax"] / (v["Tmax"] + 237.3)) - v["VP"]
        for j, expr in enumerate(exprs.values()):
            X[s, :, j] = eval(expr.replace("^", "**"), {"__builtins__": {}},
                              {**funs, **v})
    X = (X - X.mean(axis=(0, 1))) / X.std(axis=(0, 1))
    return SimpleNamespace(array=X, terms=list(exprs))


# Appends extra terms to a library.
def bind_terms(Theta, extra):
    return SimpleNamespace(array=np.concatenate([Theta.array, extra.array], axis=2),
                           terms=Theta.terms + extra.terms)


# ---- Trajectories: array Y[genotype, year, day], NaN after maturity ------
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
# Rows: all days of all years. X: library at day t. dy: y(t + 1) - y(t).
def make_design(Theta, tr, genotype):
    g = list(tr.genotypes).index(genotype)
    n = tr.nend[g]
    years = range(len(tr.years))
    X = np.vstack([Theta.array[s, :n[s], :] for s in years])
    dy = np.concatenate([np.diff(tr.Y[g, s, :n[s] + 1]) for s in years])
    year = np.repeat(np.arange(len(tr.years)), n)
    return SimpleNamespace(X=X, dy=dy, year=year, terms=Theta.terms)


# ---- Fitting methods: fit(X, dy) -> coefficients (intercept first) ---------
def fit_ols(X, dy):
    return np.linalg.lstsq(np.column_stack([np.ones(len(dy)), X]), dy,
                           rcond=None)[0]


# Lasso and elastic net, with the penalty scaling of glmnet.
def fit_glmnet(X, dy, alpha, lam):
    sx = X.std(axis=0)
    l1, l2 = lam * alpha, lam * (1 - alpha) / dy.std()
    model = ElasticNet(alpha=l1 + l2, l1_ratio=l1 / (l1 + l2), tol=1e-12,
                       max_iter=100000, precompute=True)
    model.fit((X - X.mean(axis=0)) / sx, dy - dy.mean())
    b = model.coef_ / sx
    return np.concatenate([[dy.mean() - X.mean(axis=0) @ b], b])


# Ridge: closed form, same scaling.
def fit_ridge(X, dy, lam=1e-3):
    sx = X.std(axis=0)
    Z = (X - X.mean(axis=0)) / sx
    A = Z.T @ Z / len(dy) + lam / dy.std() * np.eye(X.shape[1])
    b = np.linalg.solve(A, Z.T @ (dy - dy.mean()) / len(dy)) / sx
    return np.concatenate([[dy.mean() - X.mean(axis=0) @ b], b])


def fit_lasso(X, dy, lam=3e-4):
    return fit_glmnet(X, dy, 1, lam)


def fit_enet(X, dy, lam=6e-4):
    return fit_glmnet(X, dy, 0.5, lam)


# Best subset or backward elimination with k terms, from the Gram matrix.
def fit_subset(X, dy, k=8, method="exhaustive"):
    Xc, yc = X - X.mean(axis=0), dy - dy.mean()
    G, c = Xc.T @ Xc, Xc.T @ yc
    if method == "exhaustive":
        subsets = np.array(list(combinations(range(X.shape[1]), k)))
        blocks = G[subsets[:, :, None], subsets[:, None, :]]
        b = np.linalg.solve(blocks, c[subsets][..., None])[..., 0]
        keep = subsets[np.argmax((b * c[subsets]).sum(axis=1))]
    else:
        keep = np.arange(X.shape[1])
        while len(keep) > k:
            inv = np.linalg.inv(G[np.ix_(keep, keep)])
            increase = (inv @ c[keep]) ** 2 / np.diag(inv)
            keep = np.delete(keep, np.argmin(increase))
    b = np.zeros(X.shape[1])
    b[keep] = np.linalg.solve(G[np.ix_(keep, keep)], c[keep])
    return np.concatenate([[dy.mean() - X.mean(axis=0) @ b], b])


# Sequential thresholded least squares: refit, drop terms below the threshold.
def fit_stlsq(X, dy, threshold=0.05):
    size = X.std(axis=0, ddof=1) / dy.std(ddof=1)
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
    # "best subset (k = 8)": exhaustive search, about a second per fit
    "backward (k = 8)": lambda X, dy: fit_subset(X, dy, 8, "backward"),
    "STLSQ": fit_stlsq,
}


# ---- Simulation -------------------------------------------------------------
# Solves dy/dt = (1, theta(y, u)) coef. coef: year x (1 + terms), U: year
# x day x variable. terms(y, u) is the library; by default U is the library.
def simulate(coef, U, y0, terms=lambda y, u: u, method="LSODA"):
    D = U.shape[1]

    def rhs(t, y):
        pos = min(max(t - 0.5, 0.0), D - 1.0)   # day j is placed at time j + 0.5
        j = min(int(pos), D - 2)
        w = pos - j                             # linear interpolation
        u = (1 - w) * U[:, j] + w * U[:, j + 1]
        return coef[:, 0] + np.einsum("sk,sk->s", terms(y, u), coef[:, 1:])

    return solve_ivp(rhs, (0, D - 1), y0, method=method, t_eval=np.arange(D),
                     rtol=1e-4, atol=1e-4).y


# n genotypes spread evenly over the ranking of mean final biomass.
def spread_genotypes(tr, n=3):
    g, s = np.indices(tr.nend.shape)
    final = tr.Y[g, s, tr.nend].mean(axis=1)
    ranked = tr.genotypes[np.argsort(final)]
    return ranked[np.round(np.linspace(0, len(ranked) - 1, n)).astype(int)]


# ---- Validation: leave-one-year-out --------------------------------------
# Each year is predicted by a model fitted to the other years of the same
# genotype, starting from the mean initial biomass of those years.
# Returns P[genotype, year, day].
def cross_validate(Theta, tr, fit, genotypes=None):
    P = np.full(tr.Y.shape, np.nan)
    S = len(tr.years)
    for g in (tr.genotypes if genotypes is None else genotypes):
        i = list(tr.genotypes).index(g)
        d = make_design(Theta, tr, g)
        coef = np.array([fit(d.X[d.year != s], d.dy[d.year != s])
                         for s in range(S)])
        y0 = (tr.Y[i, :, 0].sum() - tr.Y[i, :, 0]) / (S - 1)
        P[i] = simulate(coef, Theta.array, y0)
    return P


# RMSE over whole trajectories, on a fixed day (days after sowing) and at
# maturity. Genotypes that were not predicted are ignored.
def evaluate(P, tr, fixed_day=113):
    err = P - tr.Y
    err[:, :, 0] = np.nan                        # the initial value is given
    g, s = np.indices(tr.nend.shape)

    def rmse(e):
        return np.sqrt(np.nanmean(e ** 2))

    return {"trajectory_rmse": rmse(err),
            "fixed_day_rmse": rmse(err[:, :, fixed_day - tr.das0]),
            "maturity_rmse": rmse(err[g, s, tr.nend])}


# ---- Genotype ranking -------------------------------------------------------
# Pd, Yd: genotype x year matrices of predicted and APSIM biomass on one
# day. For each year: Spearman correlation of the two rankings and overlap
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


# ---- Interpretation --------------------------------------------------------
# Standardised drivers on the days that occur in the data (year x day rows).
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
