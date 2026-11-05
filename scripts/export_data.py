"""Export the Narrabri APSIM wheat data (199 genotypes x 31 years, daily) to plain CSV files for the R workshop.
Source arrays: pool199/Narrabri.npz (derived from the public APSIM wheat simulations of Bustos-Korts et al., 2019).
Day index d=0 is day T0=15 after sowing, so days-after-sowing = T0 + d. Biomass in t/ha. Not needed to follow the workshop."""
import numpy as np, pandas as pd, gzip, sys
src, parf, out = sys.argv[1], sys.argv[2], sys.argv[3]
z = dict(np.load(src, allow_pickle=True)); T0 = int(z["T0"]); DRV = list(z["DRV"]); geno = np.asarray(z["geno"]).astype(str); years = z["years"]; tg, ty, NEND = z["tg"], z["ty"], z["NEND"]
cols = ["TT", "DL", "RAW", "RAD"]; ci = [DRV.index(c) for c in cols]
# weather: identical across genotypes within a year -> take genotype 0
wx = []
for y in range(len(years)):
    r = np.where((tg == tg.min()) & (ty == y))[0][0]; n = int(NEND[ty == y].max()); d = np.arange(n + 1)
    wx.append(pd.DataFrame(dict(year=int(years[y]), das=T0 + d, **{c: z["U"][r, :n + 1, i].astype(float) for c, i in zip(cols, ci)})))
wx = pd.concat(wx); wx["TT"] = wx.TT.round(2); wx["DL"] = wx.DL.round(3); wx["RAW"] = wx.RAW.round(4); wx["RAD"] = wx.RAD.round(2); wx.to_csv(f"{out}/narrabri_weather.csv", index=False)
# biomass and LAI up to each genotype's own maturity day
rows = []
for r in range(len(tg)):
    n = int(NEND[r]); d = np.arange(n + 1); rows.append(pd.DataFrame(dict(genotype=geno[tg[r]], year=int(years[ty[r]]), das=T0 + d, biomass=z["B"][r, :n + 1].round(5), lai=z["L"][r, :n + 1].round(4))))
bio = pd.concat(rows, ignore_index=True).sort_values(["genotype", "year", "das"]); bio.to_csv(f"{out}/narrabri_biomass.csv.gz", index=False, compression="gzip")
small = sorted(bio.genotype.unique())[::8][:25]; bio[bio.genotype.isin(small)].to_csv(f"{out}/narrabri_biomass_small.csv", index=False)
# APSIM cultivar parameters (vector parameters reduced to their first non-zero value)
P = pd.read_csv(parf); P["Genotype"] = P.Genotype.str.strip("'")
first = lambda s: [float(x) for x in str(s).strip("'").split() if float(x) != 0][0]
for c in ["transp_eff_cf", "y_extinct_coef", "y_frac_leaf", "y_rue"]: P[c] = P[c].map(first)
P = P.rename(columns={"Genotype": "genotype"}); P.round(6).to_csv(f"{out}/apsim_parameters.csv", index=False)
print("weather", wx.shape, "biomass", bio.shape, "small genotypes", len(small), "params", P.shape, "years", years[0], years[-1], "T0", T0)
