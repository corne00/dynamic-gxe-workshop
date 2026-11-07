// Compares js/builder.js with the R results of scripts/builder_reference.R.
const fs = require("fs");
const Builder = require("../js/builder.js");
const data = JSON.parse(fs.readFileSync("data/interactive/builder.json"));
const ref = JSON.parse(fs.readFileSync("scripts/builder_reference.json"));
for (const cf of ref) {
  const values = Builder.methods(1000)[cf.method].values;
  for (const r of cf.results) {
    const g = data.genotypes.indexOf(r.genotype);
    const trait = cf.trait ?? null;
    const extras = cf.extras ? Object.entries(cf.extras).map(([name, expr]) => ({name, expr})) : [];
    if (cf.states) {
      const out = Builder.runModel(data, g, {states: [].concat(cf.states), feedback: cf.feedback, drivers: trait ? [...cf.drivers, trait] : cf.drivers,
                                             degree: cf.degree, pairs: cf.pairs, method: cf.method, value: cf.value ?? null, trait, extras: []});
      console.log(cf.method.padEnd(12), [].concat(cf.states).join("+").padEnd(13), cf.feedback ? "feedback" : "        ", cf.pairs ? "pairs" : "     ",
        "R", [].concat(r.rmse).map((v) => v.toFixed(4)).join("/"), "JS", out.equations.map((e) => e.rmse.toFixed(4)).join("/"),
        "| maturity R", [].concat(r.maturity).map((v) => v.toFixed(3)).join("/"), "JS", out.equations.map((e) => e.maturityRmse.toFixed(3)).join("/"),
        "| diverged R", r.diverged, "JS", out.diverged, "|", out.ms.toFixed(0), "ms");
      continue;
    }
    const out = Builder.run(data, g, {drivers: trait ? [...cf.drivers, trait] : cf.drivers, degree: cf.degree,
                                      pairs: cf.pairs, method: cf.method, value: cf.value ?? null,
                                      target: cf.target ?? "biomass", trait, extras});
    console.log(cf.method.padEnd(12), cf.drivers.length, cf.degree, cf.pairs ? "pairs" : "     ", (cf.target ?? "biomass").padEnd(7), (cf.trait ?? "-").padEnd(7),
      "R", r.rmse.toFixed(4), "JS", out.rmse.toFixed(4),
      "| maturity R", r.maturity.toFixed(3), "JS", out.maturityRmse.toFixed(3),
      "| terms", r.nterms, out.nterms, "|", out.ms.toFixed(0), "ms");
  }
}
