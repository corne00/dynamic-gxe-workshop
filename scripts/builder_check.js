// Compares js/builder.js with the R results of scripts/builder_reference.R.
const fs = require("fs");
const Builder = require("../js/builder.js");
const data = JSON.parse(fs.readFileSync("data/interactive/builder.json"));
const ref = JSON.parse(fs.readFileSync("scripts/builder_reference.json"));
for (const cf of ref) {
  const values = Builder.methods(1000)[cf.method].values;
  for (const r of cf.results) {
    const g = data.genotypes.indexOf(r.genotype);
    const out = Builder.run(data, g, {drivers: cf.drivers, degree: cf.degree, pairs: cf.pairs,
                                      method: cf.method, value: cf.value ?? null});
    console.log(cf.method.padEnd(12), cf.drivers.length, cf.degree, cf.pairs ? "pairs" : "     ",
      "R", r.rmse.toFixed(4), "JS", out.rmse.toFixed(4),
      "| maturity R", r.maturity.toFixed(3), "JS", out.maturityRmse.toFixed(3),
      "| terms", r.nterms, out.nterms, "|", out.ms.toFixed(0), "ms");
  }
}
