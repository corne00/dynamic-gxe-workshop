// Fitting engine for the "build your own model" page. The computations are the
// ones of R/functions.R (same library, same fitting methods, leave-one-year-out
// validation), written on the Gram matrix of the design so that all years are
// fitted in milliseconds: the fit that leaves out year v uses G - G_v.
const Builder = (() => {
  // ---- Library ---------------------------------------------------------------
  function termSpecs(selected, degree, pairs) {
    const name = (x, p) => (p === 1 ? x : `${x}^${p}`);
    const specs = [], names = [];
    selected.forEach((d, i) => {
      for (let p = 1; p <= degree; p++) { specs.push([[i, p]]); names.push(name(d, p)); }
    });
    if (pairs) {
      for (let i = 0; i < selected.length; i++) {
        for (let j = i + 1; j < selected.length; j++) {
          for (let p = 1; p <= 2; p++) {
            for (let q = 1; q <= 2; q++) {
              specs.push([[i, p], [j, q]]);
              names.push(`${name(selected[i], p)}*${name(selected[j], q)}`);
            }
          }
        }
      }
    }
    return {specs, names};
  }

  // ---- Linear algebra on small dense matrices (row-major Float64Array) --------
  // Cholesky factor of a symmetric positive definite matrix; a small multiple of
  // the diagonal is added if the matrix is numerically singular.
  function cholesky(A, n) {
    let jitter = 0;
    for (let attempt = 0; attempt < 6; attempt++) {
      const L = new Float64Array(n * n);
      let ok = true;
      for (let i = 0; i < n && ok; i++) {
        for (let j = 0; j <= i; j++) {
          let s = A[i * n + j] + (i === j ? jitter * A[i * n + i] : 0);
          for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
          if (i === j) {
            if (!(s > 1e-13 * Math.abs(A[i * n + i]))) { ok = false; break; }
            L[i * n + i] = Math.sqrt(s);
          } else L[i * n + j] = s / L[j * n + j];
        }
      }
      if (ok) return L;
      jitter = jitter ? jitter * 100 : 1e-10;
    }
    return null;
  }

  function cholSolve(L, b, n) {
    const x = Float64Array.from(b);
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < i; k++) x[i] -= L[i * n + k] * x[k];
      x[i] /= L[i * n + i];
    }
    for (let i = n - 1; i >= 0; i--) {
      for (let k = i + 1; k < n; k++) x[i] -= L[k * n + i] * x[k];
      x[i] /= L[i * n + i];
    }
    return x;
  }

  function solveSPD(A, b, n) {
    const L = cholesky(A, n);
    return L ? cholSolve(L, b, n) : new Float64Array(n);
  }

  function submatrix(A, n, idx) {
    const m = idx.length, S = new Float64Array(m * m);
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) S[i * m + j] = A[idx[i] * n + idx[j]];
    return S;
  }

  function invertSPD(A, n) {
    const L = cholesky(A, n), inv = new Float64Array(n * n), e = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      e.fill(0); e[j] = 1;
      const col = cholSolve(L, e, n);
      for (let i = 0; i < n; i++) inv[i * n + j] = col[i];
    }
    return inv;
  }

  // ---- Statistics of a design from raw sums -----------------------------------
  // raw: {n, sx[K], sy, sxx[K*K], sxy[K], syy}. Returns centred quantities.
  function centred(raw, K) {
    const {n} = raw, mx = new Float64Array(K), my = raw.sy / n;
    for (let j = 0; j < K; j++) mx[j] = raw.sx[j] / n;
    const Sxx = new Float64Array(K * K), Sxy = new Float64Array(K);
    for (let j = 0; j < K; j++) {
      Sxy[j] = raw.sxy[j] - n * mx[j] * my;
      for (let k = 0; k < K; k++) Sxx[j * K + k] = raw.sxx[j * K + k] - n * mx[j] * mx[k];
    }
    return {n, K, mx, my, Sxx, Sxy, Syy: raw.syy - n * my * my};
  }

  function withIntercept(st, b) {
    let b0 = st.my;
    for (let j = 0; j < st.K; j++) b0 -= st.mx[j] * b[j];
    return {b0, b};
  }

  // ---- Fitting methods: statistics -> {b0, b} --------------------------------
  const fits = {
    OLS: (st) => withIntercept(st, solveSPD(st.Sxx, st.Sxy, st.K)),

    ridge: (st, lambda) => {
      const {K, n} = st, A = Float64Array.from(st.Sxx), sy = Math.sqrt(st.Syy / n);
      const rhs = Float64Array.from(st.Sxy, (v) => v / n);
      for (let i = 0; i < K * K; i++) A[i] /= n;
      for (let j = 0; j < K; j++) A[j * K + j] += (lambda / sy) * (st.Sxx[j * K + j] / n);
      return withIntercept(st, solveSPD(A, rhs, K));
    },

    // Coordinate descent on standardised terms, with the penalty scaling of glmnet.
    glmnet: (st, lambda, alpha, warm) => {
      const {K, n} = st, sy = Math.sqrt(st.Syy / n);
      const sx = Float64Array.from({length: K}, (_, j) => Math.sqrt(st.Sxx[j * K + j] / n));
      const l1 = lambda * alpha, l2 = (lambda * (1 - alpha)) / sy;
      const c = Float64Array.from({length: K}, (_, j) => st.Sxy[j] / (n * sx[j]));
      const G = new Float64Array(K * K);
      for (let j = 0; j < K; j++) for (let k = 0; k < K; k++) G[j * K + k] = st.Sxx[j * K + k] / (n * sx[j] * sx[k]);
      const z = warm && warm.length === K ? Float64Array.from(warm) : new Float64Array(K);
      const Gz = new Float64Array(K);               // G z, kept up to date
      for (let j = 0; j < K; j++) if (z[j] !== 0) for (let i = 0; i < K; i++) Gz[i] += G[i * K + j] * z[j];
      const update = (j) => {
        const rho = c[j] - Gz[j] + z[j];
        const next = Math.sign(rho) * Math.max(Math.abs(rho) - l1, 0) / (1 + l2);
        const delta = next - z[j];
        if (delta !== 0) {
          z[j] = next;
          for (let i = 0; i < K; i++) Gz[i] += G[i * K + j] * delta;
        }
        return Math.abs(delta);
      };
      // Full sweeps alternate with sweeps over the non-zero terms only.
      for (let outer = 0; outer < 500; outer++) {
        let change = 0;
        for (let j = 0; j < K; j++) change = Math.max(change, update(j));
        if (change < 1e-8) break;
        for (let inner = 0; inner < 2000; inner++) {
          let c2 = 0;
          for (let j = 0; j < K; j++) if (z[j] !== 0) c2 = Math.max(c2, update(j));
          if (c2 < 1e-8) break;
        }
      }
      const result = withIntercept(st, Float64Array.from(z, (v, j) => v / sx[j]));
      result.warm = z;
      return result;
    },

    // Sequential thresholded least squares.
    STLSQ: (st, threshold) => {
      const {K} = st;
      const size = Float64Array.from({length: K}, (_, j) => Math.sqrt(st.Sxx[j * K + j] / st.Syy));
      let active = Array.from({length: K}, (_, j) => j), b = new Float64Array(K);
      for (let it = 0; it < K + 1; it++) {
        const coef = solveSPD(submatrix(st.Sxx, K, active), active.map((j) => st.Sxy[j]), active.length);
        b = new Float64Array(K);
        active.forEach((j, i) => (b[j] = coef[i]));
        const keep = active.filter((j) => Math.abs(b[j]) * size[j] >= threshold);
        if (keep.length === active.length || keep.length === 0) break;
        active = keep;
      }
      return withIntercept(st, b);
    },

    // Backward elimination to k terms; the inverse is downdated at each removal.
    backward: (st, k) => {
      const {K} = st;
      let keep = Array.from({length: K}, (_, j) => j);
      let inv = invertSPD(st.Sxx, K);
      while (keep.length > k) {
        const m = keep.length;
        const coef = new Float64Array(m);
        for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) coef[i] += inv[i * m + j] * st.Sxy[keep[j]];
        let worst = 0, least = Infinity;
        for (let i = 0; i < m; i++) {
          const increase = (coef[i] * coef[i]) / inv[i * m + i];
          if (increase < least) { least = increase; worst = i; }
        }
        const next = new Float64Array((m - 1) * (m - 1));
        for (let i = 0, a = 0; i < m; i++) {
          if (i === worst) continue;
          for (let j = 0, bq = 0; j < m; j++) {
            if (j === worst) continue;
            next[a * (m - 1) + bq] = inv[i * m + j] - (inv[i * m + worst] * inv[worst * m + j]) / inv[worst * m + worst];
            bq++;
          }
          a++;
        }
        inv = next;
        keep = keep.filter((_, i) => i !== worst);
      }
      const coef = solveSPD(submatrix(st.Sxx, K, keep), keep.map((j) => st.Sxy[j]), keep.length);
      const b = new Float64Array(K);
      keep.forEach((j, i) => (b[j] = coef[i]));
      return withIntercept(st, b);
    },

    // Best subset with k terms (all subsets; limited to small searches).
    subset: (st, k) => {
      const {K} = st;
      let best = -Infinity, bestSet = null;
      const idx = [];
      const visit = (start) => {
        if (idx.length === k) {
          const c = idx.map((j) => st.Sxy[j]);
          const x = solveSPD(submatrix(st.Sxx, K, idx), c, k);
          let gain = 0;
          for (let i = 0; i < k; i++) gain += x[i] * c[i];
          if (gain > best) { best = gain; bestSet = idx.slice(); }
          return;
        }
        for (let j = start; j <= K - (k - idx.length); j++) { idx.push(j); visit(j + 1); idx.pop(); }
      };
      visit(0);
      const coef = solveSPD(submatrix(st.Sxx, K, bestSet), bestSet.map((j) => st.Sxy[j]), k);
      const b = new Float64Array(K);
      bestSet.forEach((j, i) => (b[j] = coef[i]));
      return withIntercept(st, b);
    }
  };

  function choose(n, k) {
    let r = 1;
    for (let i = 1; i <= k; i++) r = (r * (n - k + i)) / i;
    return r;
  }

  // Method table: label, name of the hyperparameter, its values, and the fit.
  const methods = (K) => {
    const kmax = Math.min(K, 20);
    const ks = Array.from({length: kmax}, (_, i) => i + 1);
    const lam = (a, b) => Array.from({length: 11}, (_, i) => 10 ** (a + ((b - a) * i) / 10));
    return {
      "OLS": {param: null, values: [null], fit: (st) => fits.OLS(st)},
      "ridge": {param: "lambda", values: lam(-5, -1), fit: (st, v) => fits.ridge(st, v)},
      "lasso": {param: "lambda", values: lam(-5, -2.5), fit: (st, v, w) => fits.glmnet(st, v, 1, w)},
      "elastic net": {param: "lambda", values: lam(-5, -2.5), fit: (st, v, w) => fits.glmnet(st, v, 0.5, w)},
      "best subset": {param: "terms", values: ks.filter((k) => choose(K, k) <= 5000),
        fit: (st, v) => fits.subset(st, v)},
      "backward": {param: "terms", values: ks, fit: (st, v) => fits.backward(st, v)},
      "STLSQ": {param: "threshold", values: [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3],
        fit: (st, v) => fits.STLSQ(st, v)}
    };
  };

  // ---- Design of one genotype ----------------------------------------------
  const cache = new Map();

  // Library rows (days 0..nend of every year), and raw sums per year.
  function prepare(data, g, selected, degree, pairs) {
    const key = [g, selected.join(","), degree, pairs].join("|");
    if (cache.has(key)) return cache.get(key);
    if (cache.size > 6) cache.delete(cache.keys().next().value);
    const {specs, names} = termSpecs(selected, degree, pairs);
    const K = specs.length, S = data.years.length;
    const col = selected.map((d) => data.drivers.indexOf(d));
    const years = [];
    const lib = (y, d) => {
      const u = col.map((c) => data.U[y][d][c]);
      return Float64Array.from(specs, (s) => s.reduce((acc, [i, p]) => acc * u[i] ** p, 1));
    };
    for (let y = 0; y < S; y++) {
      const nend = data.nend[g][y], Y = data.Y[g][y];
      const rows = Array.from({length: nend + 1}, (_, d) => lib(y, d));
      const raw = {n: nend, sx: new Float64Array(K), sy: 0, sxx: new Float64Array(K * K),
                   sxy: new Float64Array(K), syy: 0};
      for (let d = 0; d < nend; d++) {
        const x = rows[d], dy = Y[d + 1] - Y[d];
        raw.sy += dy; raw.syy += dy * dy;
        for (let j = 0; j < K; j++) {
          const xj = x[j];
          raw.sx[j] += xj; raw.sxy[j] += xj * dy;
          for (let k = 0; k <= j; k++) raw.sxx[j * K + k] += xj * x[k];
        }
      }
      for (let j = 0; j < K; j++) for (let k = j + 1; k < K; k++) raw.sxx[j * K + k] = raw.sxx[k * K + j];
      years.push({rows, raw, Y, nend});
    }
    const total = {n: 0, sx: new Float64Array(K), sy: 0, sxx: new Float64Array(K * K), sxy: new Float64Array(K), syy: 0};
    const add = (a, b, sign) => {
      a.n += sign * b.n; a.sy += sign * b.sy; a.syy += sign * b.syy;
      for (let j = 0; j < K; j++) { a.sx[j] += sign * b.sx[j]; a.sxy[j] += sign * b.sxy[j]; }
      for (let i = 0; i < K * K; i++) a.sxx[i] += sign * b.sxx[i];
    };
    years.forEach((yr) => add(total, yr.raw, 1));
    const out = {names, K, years, total, add};
    cache.set(key, out);
    return out;
  }

  // ---- Prediction: the equation integrated over a year -------------------------
  // The drivers are placed at the middle of each day and interpolated linearly, so
  // for a library of the drivers only the solution is an exact integral of the
  // piecewise linear growth rate.
  function integrate(rate, y0, nend) {
    const y = new Float64Array(nend + 1);
    y[0] = y0;
    let acc = 0.5 * rate[0];                      // from t = 0 to t = 1/2
    for (let k = 1; k <= nend; k++) {
      const mid = 0.5 * (rate[k - 1] + rate[k]);  // rate at t = k
      const piece = 0.5 * 0.5 * (rate[k - 1] + mid);   // from k - 1/2 to k
      y[k] = y0 + acc + piece;
      acc += 0.5 * (rate[k - 1] + rate[k]);        // up to t = k + 1/2
    }
    return y;
  }

  function growthRate(rows, coef, nend) {
    const K = coef.b.length, r = new Float64Array(nend + 1);
    for (let d = 0; d <= nend; d++) {
      let s = coef.b0;
      const x = rows[d];
      for (let j = 0; j < K; j++) s += x[j] * coef.b[j];
      r[d] = s;
    }
    return r;
  }

  // ---- Run: fit to all years, validate leaving one year out ------------------
  function run(data, g, config) {
    const t0 = (typeof performance !== "undefined" ? performance : Date).now();
    const design = prepare(data, g, config.drivers, config.degree, config.pairs);
    const {K, years, total} = design;
    const method = methods(K)[config.method];
    const full = method.fit(centred(total, K), config.value);
    const warm = full.warm;
    const S = years.length;
    let sq = 0, count = 0, sqEnd = 0;
    const predicted = [];
    for (let v = 0; v < S; v++) {
      const raw = {n: 0, sx: new Float64Array(K), sy: 0, sxx: new Float64Array(K * K), sxy: new Float64Array(K), syy: 0};
      design.add(raw, total, 1);
      design.add(raw, years[v].raw, -1);
      const coef = method.fit(centred(raw, K), config.value, warm);
      let y0 = 0;
      for (let u = 0; u < S; u++) if (u !== v) y0 += years[u].Y[0];
      y0 /= S - 1;
      const {rows, Y, nend} = years[v];
      const pred = integrate(growthRate(rows, coef, nend), y0, nend);
      predicted.push(Array.from(pred));
      for (let d = 1; d <= nend; d++) { sq += (pred[d] - Y[d]) ** 2; count++; }
      sqEnd += (pred[nend] - Y[nend]) ** 2;
    }
    const sd = Float64Array.from({length: K}, (_, j) => Math.sqrt(centred(total, K).Sxx[j * K + j] / (total.n - 1)));
    const sdy = Math.sqrt(centred(total, K).Syy / (total.n - 1));
    const scaled = Array.from(full.b, (v, j) => (v * sd[j]) / sdy);
    const t1 = (typeof performance !== "undefined" ? performance : Date).now();
    return {
      terms: design.names, intercept: full.b0, coef: Array.from(full.b), scaled,
      nterms: full.b.filter((v) => v !== 0).length,
      rmse: Math.sqrt(sq / count), maturityRmse: Math.sqrt(sqEnd / S),
      predicted, ms: t1 - t0
    };
  }

  return {termSpecs, methods, prepare, run, centred, fits, integrate};
})();

if (typeof module !== "undefined") module.exports = Builder;
if (typeof window !== "undefined") window.Builder = Builder;
