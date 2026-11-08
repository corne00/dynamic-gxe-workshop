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
          for (let p = 1; p < degree; p++) {
            for (let q = 1; q <= degree - p; q++) {
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
    backward: (st, k) => finish(st, backwardSet(st, k)),

    // Best subset with k terms: exact branch and bound over removals from the
    // full model. Removing term j lowers the explained sum of squares by
    // b_j^2 / inv_jj, so a branch is cut when it cannot beat the incumbent.
    subset: (st, k) => {
      const {K} = st;
      const first = backwardSet(st, k);
      const full = invertSPD(st.Sxx, K);
      const all = Array.from({length: K}, (_, j) => j);
      const explained = (keep, inv) => {
        const m = keep.length;
        let g = 0;
        for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) g += st.Sxy[keep[i]] * inv[i * m + j] * st.Sxy[keep[j]];
        return g;
      };
      const seed = first.length === K ? full : invertSPD(submatrix(st.Sxx, K, first), first.length);
      let best = explained(first, seed), bestSet = first;
      const search = (keep, inv, gain, start) => {
        const m = keep.length;
        if (m === k) return;
        const coef = new Float64Array(m);
        for (let i = 0; i < m; i++) {
          let c = 0;
          for (let j = 0; j < m; j++) c += inv[i * m + j] * st.Sxy[keep[j]];
          coef[i] = c;
        }
        const children = [];
        for (let p = start; p <= k; p++) {
          const child = gain - (coef[p] * coef[p]) / inv[p * m + p];
          if (child > best + 1e-10 * Math.abs(best)) children.push([child, p]);
        }
        children.sort((x, y) => y[0] - x[0]);
        for (const [child, p] of children) {
          if (child <= best + 1e-10 * Math.abs(best)) continue;
          const rest = keep.filter((_, i) => i !== p);
          if (rest.length === k) { best = child; bestSet = rest; continue; }
          search(rest, downdate(inv, m, p), child, p);
        }
      };
      search(all, full, explained(all, full), 0);
      return finish(st, bestSet);
    }
  };

  // Inverse of a symmetric matrix after deleting row and column p.
  function downdate(inv, m, p) {
    const next = new Float64Array((m - 1) * (m - 1));
    const d = inv[p * m + p];
    for (let i = 0, a = 0; i < m; i++) {
      if (i === p) continue;
      const f = inv[i * m + p] / d;
      for (let j = 0, c = 0; j < m; j++) {
        if (j === p) continue;
        next[a * (m - 1) + c] = inv[i * m + j] - f * inv[p * m + j];
        c++;
      }
      a++;
    }
    return next;
  }

  function backwardSet(st, k) {
    const {K} = st;
    let keep = Array.from({length: K}, (_, j) => j);
    let inv = invertSPD(st.Sxx, K);
    while (keep.length > k) {
      const m = keep.length;
      let worst = 0, least = Infinity;
      for (let i = 0; i < m; i++) {
        let c = 0;
        for (let j = 0; j < m; j++) c += inv[i * m + j] * st.Sxy[keep[j]];
        const increase = (c * c) / inv[i * m + i];
        if (increase < least) { least = increase; worst = i; }
      }
      inv = downdate(inv, m, worst);
      keep = keep.filter((_, i) => i !== worst);
    }
    return keep;
  }

  function finish(st, keep) {
    const {K} = st;
    const coef = solveSPD(submatrix(st.Sxx, K, keep), keep.map((j) => st.Sxy[j]), keep.length);
    const b = new Float64Array(K);
    keep.forEach((j, i) => (b[j] = coef[i]));
    return withIntercept(st, b);
  }

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

  // ---- Extra terms: expressions of the raw weather -----------------------------
  // A small expression language (no code is evaluated): numbers, the variables of
  // the weather table, + - * / ^, parentheses, and the functions exp, log, sqrt,
  // abs, sin, cos, tan, tanh, min, max (elementwise), roll(x, k) (mean of the last
  // k days), cum(x) (cumulative sum) and lag(x, k).
  const functions = {
    exp: Math.exp, log: Math.log, sqrt: Math.sqrt, abs: Math.abs,
    sin: Math.sin, cos: Math.cos, tan: Math.tan, tanh: Math.tanh
  };

  function tokenize(src) {
    const tokens = [], re = /\s*(?:(\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+)|([A-Za-z_][A-Za-z0-9_.]*)|(.))/gy;
    let m;
    while (re.lastIndex < src.length && (m = re.exec(src))) {
      if (m[1] !== undefined) tokens.push({type: "num", value: +m[1]});
      else if (m[2] !== undefined) tokens.push({type: "id", value: m[2]});
      else if (m[3] !== undefined && m[3].trim()) tokens.push({type: "op", value: m[3]});
    }
    return tokens;
  }

  function parseExpression(src) {
    const t = tokenize(src);
    let i = 0;
    const peek = () => t[i], take = () => t[i++];
    const expect = (v) => {
      if (!peek() || peek().value !== v) throw new Error(`expected "${v}"`);
      take();
    };
    function atom() {
      const tok = take();
      if (!tok) throw new Error("unexpected end of expression");
      if (tok.type === "num") return {num: tok.value};
      if (tok.type === "id") {
        if (peek() && peek().value === "(") {
          take();
          const args = [];
          if (peek() && peek().value !== ")") {
            args.push(sum());
            while (peek() && peek().value === ",") { take(); args.push(sum()); }
          }
          expect(")");
          return {call: tok.value, args};
        }
        return {id: tok.value};
      }
      if (tok.value === "(") { const e = sum(); expect(")"); return e; }
      throw new Error(`unexpected "${tok.value}"`);
    }
    function power() {
      const base = atom();
      if (peek() && peek().value === "^") { take(); return {op: "^", a: base, b: unary()}; }
      return base;
    }
    function unary() {
      if (peek() && peek().value === "-") { take(); return {op: "neg", a: unary()}; }
      if (peek() && peek().value === "+") { take(); return unary(); }
      return power();
    }
    function product() {
      let e = unary();
      while (peek() && (peek().value === "*" || peek().value === "/")) e = {op: take().value, a: e, b: unary()};
      return e;
    }
    function sum() {
      let e = product();
      while (peek() && (peek().value === "+" || peek().value === "-")) e = {op: take().value, a: e, b: product()};
      return e;
    }
    const ast = sum();
    if (i < t.length) throw new Error(`unexpected "${t[i].value}"`);
    return ast;
  }

  // Evaluates an expression on one year; vars maps names to Float64Array(D).
  function evaluate(ast, vars, D) {
    const vec = (v) => (typeof v === "number" ? new Float64Array(D).fill(v) : v);
    const map2 = (a, b, f) => {
      const x = vec(a), y = vec(b), out = new Float64Array(D);
      for (let d = 0; d < D; d++) out[d] = f(x[d], y[d]);
      return out;
    };
    const scalar = (node, what) => {
      if (node.num === undefined) throw new Error(`${what} must be a number`);
      return node.num;
    };
    const go = (n) => {
      if (n.num !== undefined) return n.num;
      if (n.id !== undefined) {
        if (n.id === "pi") return Math.PI;
        if (!(n.id in vars)) throw new Error(`unknown variable "${n.id}"`);
        return vars[n.id];
      }
      if (n.op === "neg") return map2(go(n.a), 0, (x) => -x);
      if (n.op) {
        const f = {"+": (x, y) => x + y, "-": (x, y) => x - y, "*": (x, y) => x * y,
                   "/": (x, y) => x / y, "^": (x, y) => x ** y}[n.op];
        return map2(go(n.a), go(n.b), f);
      }
      const f = n.call;
      if (functions[f]) {
        if (n.args.length !== 1) throw new Error(`${f} takes one argument`);
        const a = go(n.args[0]);
        return typeof a === "number" ? functions[f](a) : Float64Array.from(a, functions[f]);
      }
      if (f === "min" || f === "max") {
        if (n.args.length < 2) throw new Error(`${f} takes at least two arguments`);
        const pick = f === "min" ? Math.min : Math.max;
        return n.args.map(go).reduce((a, b) => map2(a, b, pick));
      }
      if (f === "roll" || f === "lag") {
        if (n.args.length !== 2) throw new Error(`${f} takes two arguments`);
        const x = vec(go(n.args[0])), k = scalar(n.args[1], "the second argument of " + f);
        const out = new Float64Array(D);
        if (f === "lag") for (let d = 0; d < D; d++) out[d] = x[Math.max(d - k, 0)];
        else {
          const run = new Float64Array(D + 1);
          for (let d = 0; d < D; d++) run[d + 1] = run[d] + x[d];
          for (let d = 0; d < D; d++) { const from = Math.max(0, d - k + 1); out[d] = (run[d + 1] - run[from]) / (d + 1 - from); }
        }
        return out;
      }
      if (f === "cum") {
        if (n.args.length !== 1) throw new Error("cum takes one argument");
        const x = vec(go(n.args[0])), out = new Float64Array(D);
        let a = 0;
        for (let d = 0; d < D; d++) out[d] = a += x[d];
        return out;
      }
      throw new Error(`unknown function "${f}"`);
    };
    return vec(go(ast));
  }

  // One line per term, "name = expression" or just "expression".
  const dummy = Object.fromEntries(["Tmin", "Tmax", "Rain", "RAD", "VP", "DL", "TT", "RAW", "Tmean", "VPD", "das"]
    .map((v) => [v, new Float64Array([1, 2, 3])]));
  function parseExtras(text, traits = []) {
    const terms = [], errors = [];
    const known = Object.assign({}, dummy);
    traits.forEach((t) => { known[t] = new Float64Array([1, 2, 3]); });
    text.split("\n").forEach((line, i) => {
      const raw = line.trim();
      if (!raw || raw.startsWith("#")) return;
      const eq = raw.match(/^([A-Za-z_][A-Za-z0-9_.]*)\s*=(?!=)\s*(.+)$/);
      const [name, expr] = eq ? [eq[1], eq[2]] : [raw.replace(/\s+/g, ""), raw];
      try { evaluate(parseExpression(expr), known, 3); terms.push({name, expr}); }
      catch (e) { errors.push({line: i + 1, text: raw, message: e.message}); }
    });
    return {terms, errors};
  }

  // ---- Measured states used like drivers --------------------------------------
  // A record is the value at the end of a day, so the value of a day is the mean
  // of the records at its start and end; later days repeat the last value. The
  // standardised series are scaled like the drivers (over all years and days).
  const traitCache = new Map();
  function traitSeries(data, g, name) {
    const key = `${g}|${name}`;
    if (traitCache.has(key)) return traitCache.get(key);
    const source = name === "lai" ? data.L[g] : data.Y[g];
    const D = data.U[0].length;
    const raw = source.map((rec) => {
      const day = new Float64Array(D);
      for (let d = 0; d < D; d++) {
        day[d] = d === 0 ? rec[0] : d < rec.length ? (rec[d - 1] + rec[d]) / 2 : day[d - 1];
      }
      return day;
    });
    let n = 0, sum = 0, sq = 0;
    for (const row of raw) for (const x of row) { n++; sum += x; sq += x * x; }
    const mean = sum / n, sd = Math.sqrt(sq / n - mean * mean);
    const out = {raw, std: raw.map((row) => Float64Array.from(row, (x) => (x - mean) / sd))};
    traitCache.set(key, out);
    return out;
  }

  // Standardised columns [year][day] for the extra terms (like the drivers).
  const extraCache = new Map();
  function extraColumns(data, extras, g = 0, trait = null) {
    const key = extras.map((e) => `${e.name}=${e.expr}`).join(";") + `|${g}|${trait}`;
    if (extraCache.has(key)) return extraCache.get(key);
    const S = data.years.length, D = data.U[0].length;
    const series = trait ? traitSeries(data, g, trait) : null;
    const raw = (y) => {
      const v = {};
      data.drivers.forEach((name, k) => {
        v[name] = Float64Array.from({length: D}, (_, d) => data.U[y][d][k] * data.scale[k] + data.center[k]);
      });
      v.Tmean = Float64Array.from({length: D}, (_, d) => (v.Tmin[d] + v.Tmax[d]) / 2);
      v.VPD = Float64Array.from({length: D}, (_, d) => 6.108 * Math.exp((17.27 * v.Tmax[d]) / (v.Tmax[d] + 237.3)) - v.VP[d]);
      v.das = Float64Array.from({length: D}, (_, d) => data.das0 + d);
      if (series) v[trait] = series.raw[y];
      return v;
    };
    const vars = Array.from({length: S}, (_, y) => raw(y));
    const columns = extras.map((e) => {
      const ast = parseExpression(e.expr);
      const values = vars.map((v) => evaluate(ast, v, D));
      let n = 0, sum = 0, sq = 0;
      for (const row of values) for (const x of row) {
        if (!Number.isFinite(x)) throw new Error(`"${e.name}" is not finite on some days (for example a logarithm of a negative number)`);
        n++; sum += x; sq += x * x;
      }
      const mean = sum / n, sd = Math.sqrt(Math.max(sq / n - mean * mean, 0));
      if (!(sd > 1e-12)) throw new Error(`"${e.name}" is constant`);
      return {name: e.name, values: values.map((row) => Float64Array.from(row, (x) => (x - mean) / sd))};
    });
    if (extraCache.size > 20) extraCache.delete(extraCache.keys().next().value);
    extraCache.set(key, columns);
    return columns;
  }

  // ---- Design of one genotype ----------------------------------------------
  const cache = new Map();

  // Library rows (days 0..nend of every year), and raw sums per year.
  function prepare(data, g, selected, degree, pairs, extras = [], target = "biomass", trait = null) {
    const key = [g, selected.join(","), degree, pairs, extras.map((e) => `${e.name}=${e.expr}`).join(";"), target, trait].join("|");
    if (cache.has(key)) return cache.get(key);
    if (cache.size > 6) cache.delete(cache.keys().next().value);
    const base = termSpecs(selected, degree, pairs), specs = base.specs;
    const extra = extraColumns(data, extras, g, trait);
    const series = trait ? traitSeries(data, g, trait) : null;
    const names = base.names.concat(extra.map((e) => e.name));
    const K = specs.length + extra.length, S = data.years.length;
    const col = selected.map((d) => (d === trait ? -1 : data.drivers.indexOf(d)));
    const targetSeries = target === "lai" ? data.L[g] : data.Y[g];
    const years = [];
    const lib = (y, d) => {
      const u = col.map((c) => (c < 0 ? series.std[y][d] : data.U[y][d][c]));
      const row = new Float64Array(K);
      specs.forEach((s, j) => { row[j] = s.reduce((acc, [i, p]) => acc * u[i] ** p, 1); });
      extra.forEach((e, j) => { row[specs.length + j] = e.values[y][d]; });
      return row;
    };
    for (let y = 0; y < S; y++) {
      const nend = data.nend[g][y], Y = targetSeries[y];
      const rows = Array.from({length: Math.min(nend + 2, data.U[y].length)}, (_, d) => lib(y, d));
      const raw = {n: 0, sx: new Float64Array(K), sy: 0, sxx: new Float64Array(K * K),
                   sxy: new Float64Array(K), syy: 0};
      raw.n = nend;
      const x = new Float64Array(K);
      for (let d = 1; d <= nend; d++) {
        const dy = Y[d] - Y[d - 1];                 // growth rate of day d, one-day difference
        for (let j = 0; j < K; j++) x[j] = rows[d][j];
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
  // Biomass is recorded at the end of a day, so day d (counted from 0) drives the
  // interval ending at record d. Its value is placed at time d - 1/2 (time 0 is the
  // first record) and interpolated linearly; for a library of the drivers only the
  // solution is the exact integral of this piecewise linear growth rate.
  function integrate(rate, y0, nend) {
    const r = (k) => rate[Math.min(k, rate.length - 1)];
    const y = new Float64Array(nend + 1);
    y[0] = y0;
    let acc = 0.5 * 0.5 * (0.5 * (r(0) + r(1)) + r(1));      // from t = 0 to t = 1/2
    for (let k = 1; k <= nend; k++) {
      const mid = 0.5 * (r(k) + r(k + 1));                    // rate at t = k
      y[k] = y0 + acc + 0.5 * 0.5 * (r(k) + mid);              // plus k - 1/2 to k
      acc += 0.5 * (r(k) + r(k + 1));                          // up to t = k + 1/2
    }
    return y;
  }

  function growthRate(rows, coef, nend) {
    const K = coef.b.length, r = new Float64Array(rows.length);
    for (let d = 0; d < rows.length; d++) {
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
    const design = prepare(data, g, config.drivers, config.degree, config.pairs, config.extras || [], config.target || "biomass", config.trait || null);
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

  // ---- State feedback and coupled models ---------------------------------------
  // The modelled states enter the library through standardised terms (centre and
  // standard deviation of the genotype's records): powers of each state, and with
  // products also products of states and products of a state with powers of
  // 2 of each driver (counterpart of state_terms in R/functions.R). Every state has
  // its own equation, fitted with the same method; the system is integrated with
  // fixed-step RK4 (a quarter of a day), and a held-out year whose solution exceeds
  // ten times the largest observed value is reported as diverged.
  const stateSeries = (data, g, name) => (name === "lai" ? data.L[g] : data.Y[g]);

  function stateScale(series) {
    let n = 0, sum = 0, sq = 0;
    for (const rec of series) for (const x of rec) { n++; sum += x; sq += x * x; }
    const center = sum / n;
    return {center, scale: Math.sqrt(sq / n - center * center)};
  }

  // Terms {a, p, b, q, c}: z_a^p times z_b^q (b >= 0), a library column c (c >= 0) or 1.
  function stateSpecs(states, base, degree, pairs) {
    const name = (x, p) => (p === 1 ? x : `${x}^${p}`);
    const inter = [];                     // library columns that interact with a state
    base.specs.forEach((sp, j) => { if (sp.length === 1 && sp[0][1] <= 2) inter.push(j); });
    const specs = [], names = [];
    states.forEach((a, i) => {
      for (let p = 1; p <= degree; p++) { specs.push({a: i, p, b: -1, q: 0, c: -1}); names.push(name(a, p)); }
    });
    if (pairs) {
      for (let i = 0; i < states.length; i++) for (let j = i + 1; j < states.length; j++) {
        for (let p = 1; p < degree; p++) for (let q = 1; q <= degree - p; q++) {
          specs.push({a: i, p, b: j, q, c: -1}); names.push(`${name(states[i], p)}*${name(states[j], q)}`);
        }
      }
      states.forEach((a, i) => {
        for (const c of inter) for (let p = 1; p <= degree - base.specs[c][0][1]; p++) {
          specs.push({a: i, p, b: -1, q: 0, c}); names.push(`${name(a, p)}*${base.names[c]}`);
        }
      });
    }
    return {specs, names};
  }

  // State terms for the states z (standardised) and the library columns x.
  function stateRow(specs, z, x, out, offset) {
    for (let s = 0; s < specs.length; s++) {
      const t = specs[s];
      let v = z[t.a] ** t.p;
      if (t.b >= 0) v *= z[t.b] ** t.q;
      else if (t.c >= 0) v *= x[t.c];
      out[offset + s] = v;
    }
  }

  const stateCache = new Map();
  function prepareStates(data, g, config) {
    const {states, drivers, degree, pairs, feedback, trait = null} = config, extras = config.extras || [];
    const key = [g, states.join("+"), drivers.join(","), degree, pairs, feedback, trait,
                 extras.map((e) => `${e.name}=${e.expr}`).join(";")].join("|");
    if (stateCache.has(key)) return stateCache.get(key);
    if (stateCache.size > 4) stateCache.delete(stateCache.keys().next().value);
    const base = prepare(data, g, drivers, degree, pairs, extras, states[0], trait);
    const Kw = base.K, S = data.years.length, m = states.length;
    const series = states.map((n) => stateSeries(data, g, n));
    const scales = series.map(stateScale);
    const ss = feedback ? stateSpecs(states, termSpecs(drivers, degree, pairs), degree, pairs) : {specs: [], names: []};
    const Kt = Kw + ss.specs.length, names = base.names.concat(ss.names);
    const years = [], z = new Float64Array(m), x = new Float64Array(Kt);
    const total = {n: 0, sx: new Float64Array(Kt), sxx: new Float64Array(Kt * Kt),
                   sy: new Float64Array(m), sxy: Array.from({length: m}, () => new Float64Array(Kt)), syy: new Float64Array(m)};
    for (let y = 0; y < S; y++) {
      const {rows} = base.years[y], nend = data.nend[g][y];
      const raw = {n: nend, sx: new Float64Array(Kt), sxx: new Float64Array(Kt * Kt),
                   sy: new Float64Array(m), sxy: Array.from({length: m}, () => new Float64Array(Kt)), syy: new Float64Array(m)};
      for (let d = 1; d <= nend; d++) {
        for (let j = 0; j < Kw; j++) x[j] = rows[d][j];
        if (feedback) {
          for (let k = 0; k < m; k++) z[k] = ((series[k][y][d - 1] + series[k][y][d]) / 2 - scales[k].center) / scales[k].scale;
          stateRow(ss.specs, z, x, x, Kw);
        }
        for (let k = 0; k < m; k++) {
          const dy = series[k][y][d] - series[k][y][d - 1];
          raw.sy[k] += dy; raw.syy[k] += dy * dy;
          for (let j = 0; j < Kt; j++) raw.sxy[k][j] += x[j] * dy;
        }
        for (let j = 0; j < Kt; j++) {
          raw.sx[j] += x[j];
          for (let i = 0; i <= j; i++) raw.sxx[j * Kt + i] += x[j] * x[i];
        }
      }
      for (let j = 0; j < Kt; j++) for (let i = j + 1; i < Kt; i++) raw.sxx[j * Kt + i] = raw.sxx[i * Kt + j];
      years.push({rows, raw, nend, Y: series.map((sr) => sr[y])});
    }
    const add = (a, b, sign) => {
      a.n += sign * b.n;
      for (let j = 0; j < Kt; j++) a.sx[j] += sign * b.sx[j];
      for (let i = 0; i < Kt * Kt; i++) a.sxx[i] += sign * b.sxx[i];
      for (let k = 0; k < m; k++) {
        a.sy[k] += sign * b.sy[k]; a.syy[k] += sign * b.syy[k];
        for (let j = 0; j < Kt; j++) a.sxy[k][j] += sign * b.sxy[k][j];
      }
    };
    years.forEach((yr) => add(total, yr.raw, 1));
    const limits = series.map((sr) => 10 * Math.max(...sr.map((rec) => Math.max(...rec.map(Math.abs)))));
    const out = {names, Kw, Kt, years, total, add, specs: ss.specs, scales, limits};
    stateCache.set(key, out);
    return out;
  }

  // The coupled system over one year; null if the solution blows up.
  function integrateStates(rows, coefs, design, y0, nend, feedback) {
    const {Kw, specs, scales, limits} = design, m = coefs.length, D = rows.length;
    const u = new Float64Array(Kw), full = new Float64Array(Kw + specs.length), z = new Float64Array(m);
    const rate = (t, y, out) => {
      const pos = Math.min(Math.max(t + 0.5, 0), D - 1), j = Math.min(Math.floor(pos), D - 2), w = pos - j;
      const r0 = rows[j], r1 = rows[j + 1];
      for (let c = 0; c < Kw; c++) full[c] = u[c] = (1 - w) * r0[c] + w * r1[c];
      if (feedback) {
        for (let k = 0; k < m; k++) z[k] = (y[k] - scales[k].center) / scales[k].scale;
        stateRow(specs, z, u, full, Kw);
      }
      for (let k = 0; k < m; k++) {
        let s = coefs[k].b0; const b = coefs[k].b;
        for (let c = 0; c < full.length; c++) s += b[c] * full[c];
        out[k] = s;
      }
    };
    const steps = 4, h = 1 / steps;
    const pred = Array.from({length: m}, () => new Float64Array(nend + 1));
    let y = Float64Array.from(y0);
    const k1 = new Float64Array(m), k2 = new Float64Array(m), k3 = new Float64Array(m), k4 = new Float64Array(m), tmp = new Float64Array(m);
    for (let k = 0; k < m; k++) pred[k][0] = y[k];
    for (let d = 1; d <= nend; d++) {
      for (let s = 0; s < steps; s++) {
        const t = d - 1 + s * h;
        rate(t, y, k1);
        for (let k = 0; k < m; k++) tmp[k] = y[k] + 0.5 * h * k1[k];
        rate(t + 0.5 * h, tmp, k2);
        for (let k = 0; k < m; k++) tmp[k] = y[k] + 0.5 * h * k2[k];
        rate(t + 0.5 * h, tmp, k3);
        for (let k = 0; k < m; k++) tmp[k] = y[k] + h * k3[k];
        rate(t + h, tmp, k4);
        for (let k = 0; k < m; k++) y[k] += (h / 6) * (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]);
        for (let k = 0; k < m; k++) if (!(Math.abs(y[k]) <= limits[k])) return null;
      }
      for (let k = 0; k < m; k++) pred[k][d] = y[k];
    }
    return pred;
  }

  // Fit to all years, validate leaving one year out, for one or two coupled states.
  // config: states (["biomass"], ["lai"] or both), feedback, drivers (with the
  // measured trait, if any), trait, degree, pairs, extras, method, value.
  function runModel(data, g, config) {
    const states = config.states, feedback = !!config.feedback;
    if (states.length === 1 && !feedback) {
      const r = run(data, g, {...config, target: states[0]});
      return {equations: [{state: states[0], ...r}], diverged: 0, ms: r.ms, scales: []};
    }
    const now = () => (typeof performance !== "undefined" ? performance : Date).now();
    const t0 = now();
    const design = prepareStates(data, g, {...config, states, feedback});
    const {Kt, years, total} = design, m = states.length, S = years.length;
    const method = methods(Kt)[config.method];
    const view = (raw, k) => centred({n: raw.n, sx: raw.sx, sxx: raw.sxx, sy: raw.sy[k], sxy: raw.sxy[k], syy: raw.syy[k]}, Kt);
    const full = states.map((_, k) => method.fit(view(total, k), config.value));
    const sq = new Float64Array(m), sqEnd = new Float64Array(m), counts = new Float64Array(m);
    let nEnd = 0, diverged = 0;
    const predicted = states.map(() => []);
    for (let v = 0; v < S; v++) {
      const raw = {n: 0, sx: new Float64Array(Kt), sxx: new Float64Array(Kt * Kt), sy: new Float64Array(m),
                   sxy: Array.from({length: m}, () => new Float64Array(Kt)), syy: new Float64Array(m)};
      design.add(raw, total, 1); design.add(raw, years[v].raw, -1);
      const coefs = states.map((_, k) => method.fit(view(raw, k), config.value, full[k].warm));
      const y0 = states.map((_, k) => {
        let s = 0;
        for (let u = 0; u < S; u++) if (u !== v) s += years[u].Y[k][0];
        return s / (S - 1);
      });
      const {rows, Y, nend} = years[v];
      const pred = integrateStates(rows, coefs, design, y0, nend, feedback);
      if (!pred) {
        diverged++;
        states.forEach((_, k) => predicted[k].push(null));
        continue;
      }
      nEnd++;
      states.forEach((_, k) => {
        predicted[k].push(Array.from(pred[k]));
        for (let d = 1; d <= nend; d++) { sq[k] += (pred[k][d] - Y[k][d]) ** 2; counts[k]++; }
        sqEnd[k] += (pred[k][nend] - Y[k][nend]) ** 2;
      });
    }
    const equations = states.map((state, k) => {
      const st = view(total, k);
      const sd = Float64Array.from({length: Kt}, (_, j) => Math.sqrt(st.Sxx[j * Kt + j] / (total.n - 1)));
      const sdy = Math.sqrt(st.Syy / (total.n - 1));
      return {state, terms: design.names, intercept: full[k].b0, coef: Array.from(full[k].b),
              scaled: Array.from(full[k].b, (b, j) => (b * sd[j]) / sdy),
              nterms: full[k].b.filter((b) => b !== 0).length,
              rmse: Math.sqrt(sq[k] / counts[k]), maturityRmse: Math.sqrt(sqEnd[k] / nEnd),
              predicted: predicted[k]};
    });
    return {equations, diverged, ms: now() - t0,
            scales: feedback ? states.map((state, k) => ({state, ...design.scales[k]})) : []};
  }

  // ---- Views of the design (steps 2 and 3 of the build page) -------------------
  // Correlations between the library columns, and between each column and the
  // growth rate of each modelled state, over all rows of the genotype's design.
  function correlations(data, g, config) {
    const states = config.states, feedback = !!config.feedback;
    let names, K, Sxx, n;
    const Sxy = [], Syy = [];
    if (states.length === 1 && !feedback) {
      const d = prepare(data, g, config.drivers, config.degree, config.pairs, config.extras || [], states[0], config.trait || null);
      const c = centred(d.total, d.K);
      names = d.names; K = d.K; Sxx = c.Sxx; n = c.n;
      Sxy.push(c.Sxy); Syy.push(c.Syy);
    } else {
      const d = prepareStates(data, g, {...config, states, feedback});
      names = d.names; K = d.Kt; n = d.total.n;
      states.forEach((_, k) => {
        const c = centred({n, sx: d.total.sx, sxx: d.total.sxx, sy: d.total.sy[k], sxy: d.total.sxy[k], syy: d.total.syy[k]}, K);
        Sxx = c.Sxx; Sxy.push(c.Sxy); Syy.push(c.Syy);
      });
    }
    const sd = Float64Array.from({length: K}, (_, j) => Math.sqrt(Math.max(Sxx[j * K + j], 0)));
    const corr = new Float64Array(K * K);
    for (let i = 0; i < K; i++) for (let j = 0; j < K; j++) {
      const den = sd[i] * sd[j];
      corr[i * K + j] = den > 1e-12 ? Sxx[i * K + j] / den : 0;
    }
    const withY = Sxy.map((sxy, k) => Float64Array.from(sxy, (v, j) => {
      const den = sd[j] * Math.sqrt(Math.max(Syy[k], 0));
      return den > 1e-12 ? v / den : 0;
    }));
    return {names, K, n, corr, withY};
  }

  // Rows of the design matrix of one year: days from .. from + count - 1 (counted
  // as records, the first record of a year has no row), with the target.
  function designRows(data, g, config, y, from, count) {
    const d = prepare(data, g, config.drivers, config.degree, config.pairs, config.extras || [], config.states[0], config.trait || null);
    const yr = d.years[y], rows = [];
    for (let k = from; k < Math.min(from + count, yr.nend + 1); k++) {
      rows.push({day: k, x: Array.from(yr.rows[k]), dy: yr.Y[k] - yr.Y[k - 1], y0: yr.Y[k - 1], y1: yr.Y[k]});
    }
    return {names: d.names, rows, nend: yr.nend};
  }

  // One library column against the growth rate, all rows of all years.
  function designColumn(data, g, config, j) {
    const d = prepare(data, g, config.drivers, config.degree, config.pairs, config.extras || [], config.states[0], config.trait || null);
    const out = [];
    d.years.forEach((yr, y) => {
      for (let k = 1; k <= yr.nend; k++) out.push({year: y, day: k, x: yr.rows[k][j], dy: yr.Y[k] - yr.Y[k - 1]});
    });
    return out;
  }

  return {termSpecs, methods, prepare, run, runModel, stateSpecs, prepareStates, centred, fits, integrate, parseExtras, parseExpression, evaluate, extraColumns, traitSeries, correlations, designRows, designColumn};
})();

if (typeof module !== "undefined") module.exports = Builder;
if (typeof window !== "undefined") window.Builder = Builder;
