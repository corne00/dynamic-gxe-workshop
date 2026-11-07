// The fitted equation as LaTeX: dy/dt = c0 + c1 term1 + ..., typeset with MathJax.
const Equation = (() => {
  const pretty = {lai: "LAI", biomass: "B"};
  const word = (name) => `\\mathrm{${(pretty[name] || name).replace(/_/g, "\\_")}}`;
  // "TT^2*RAW" -> \mathrm{TT}^{2}\,\mathrm{RAW}
  function term(name) {
    if (/[^A-Za-z0-9_.^*]/.test(name)) return `\\mathrm{${name.replace(/[_#%&$]/g, "\\$&")}}`;
    return name.split("*").map((f) => {
      const [base, power] = f.split("^");
      return power ? `${word(base)}^{${power}}` : word(base);
    }).join("\\,");
  }
  function number(v) {
    const a = Math.abs(v);
    if (a >= 0.01 && a < 1000) return `${+a.toPrecision(3)}`;
    const [m, e] = a.toExponential(2).split("e");
    return `${+m}\\times 10^{${+e}}`;
  }
  // Terms with a non-zero coefficient, `per` terms to a line.
  function latex(terms, intercept, coef, lhs = "\\frac{dy}{dt}", per = 3) {
    const parts = [];
    if (intercept !== 0) parts.push({sign: intercept < 0 ? "-" : "+", text: number(intercept)});
    coef.forEach((c, i) => { if (c !== 0) parts.push({sign: c < 0 ? "-" : "+", text: `${number(c)}\\,${term(terms[i])}`}); });
    if (!parts.length) return `${lhs} = 0`;
    const lines = [];
    for (let i = 0; i < parts.length; i += per) {
      lines.push(parts.slice(i, i + per).map((p, k) => {
        const first = i === 0 && k === 0;
        return first ? (p.sign === "-" ? `-${p.text}` : p.text) : ` ${p.sign} ${p.text}`;
      }).join(""));
    }
    return `\\begin{aligned} ${lhs} &= ${lines.join(" \\\\ &\\quad ")} \\end{aligned}`;
  }
  // A display equation element; typeset once MathJax is available.
  function element(source) {
    const el = document.createElement("div");
    el.className = "dgxe-equation";
    el.textContent = `\\[${source}\\]`;
    const typeset = (tries) => {
      const mj = window.MathJax;
      if (mj && mj.typesetPromise) mj.typesetPromise([el]).catch(() => {});
      else if (tries > 0) setTimeout(() => typeset(tries - 1), 200);
    };
    typeset(25);
    return el;
  }
  return {latex, element, term, number};
})();
if (typeof module !== "undefined") module.exports = Equation;
if (typeof window !== "undefined") window.Equation = Equation;
