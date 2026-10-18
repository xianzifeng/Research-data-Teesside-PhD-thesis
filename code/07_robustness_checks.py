"""
07_robustness_checks.py  (thesis Appendix C.5; Section 5.5.1; Tables 5-5, 5-8 to 5-12; Appendix D.3, D.4)
  - leave-one-year-out re-estimation (Table 5-9)
  - alternative measurement specification: equal-weight vs loading-weighted composites (Table 5-10)
  - placebo test with permuted FA (1,000 permutations)
  - leave-one-unit-out (Table 5-12)
  - outlier treatments: winsorising at 1 %, dropping Cook's D flags (Table D-4)
  - regional subgroups (Table D-3)
  - wild-cluster bootstrap by business unit for the few-clusters problem (comment C476)
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

from utils import OUT, BLOCKS, SECOND_ORDER, banner, standardized_loadings

df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv"
banner("07 Robustness checks", df)
cl = lambda d: {"cov_type": "cluster", "cov_kwds": {"groups": d["employee_id"]}}

def fit(d, f="zDC ~ zFA + zRR"):
    r = smf.ols(f, data=d).fit(**cl(d))
    return r.params["zFA"], r.params["zRR"]

base = fit(df)
print(f"Baseline composite model: FA->DC = {base[0]:.3f}, RR->DC = {base[1]:.3f}")

# 1. Leave-one-year-out (Table 5-9)
print("\nLeave-one-year-out:")
for w in sorted(df.wave.unique()):
    b = fit(df[df.wave != w]); print(f"   excluding wave {w} ({2019 + w}): FA->DC = {b[0]:.3f}, RR->DC = {b[1]:.3f}")

# 2. Alternative measurement specification (Table 5-10): loading-weighted composites
d2 = df.copy()
for blk, cols in BLOCKS.items():
    lam = standardized_loadings(df[cols]); d2[blk + "_w"] = (df[cols] * lam).sum(axis=1) / lam.sum()
d2["RR_w"] = d2[[b + "_w" for b in SECOND_ORDER["RR"]]].mean(axis=1)
d2["DC_w"] = d2[[b + "_w" for b in SECOND_ORDER["DC"]]].mean(axis=1)
for c in ["FA_w", "RR_w", "DC_w"]:
    d2["z" + c] = (d2[c] - d2[c].mean()) / d2[c].std()
a = smf.ols("zRR_w ~ zFA_w", data=d2).fit(**cl(d2)).params["zFA_w"]
o = smf.ols("zDC_w ~ zFA_w + zRR_w", data=d2).fit(**cl(d2)).params
print(f"\nLoading-weighted composites: a = {a:.3f}, b = {o['zRR_w']:.3f}, c' = {o['zFA_w']:.3f}, indirect = {a*o['zRR_w']:.3f}")

# 3. Placebo (permuted FA), 1,000 permutations
rng = np.random.default_rng(42)
pl = []
for _ in range(1000):
    d = df.assign(FA_placebo=rng.permutation(df.zFA.values))
    pl.append(sm.OLS(d.zDC, sm.add_constant(d[["FA_placebo", "zRR"]])).fit().params["FA_placebo"])
print(f"Placebo FA->DC 95% band: [{np.percentile(pl, 2.5):.3f}, {np.percentile(pl, 97.5):.3f}]  (observed {base[0]:.3f})")

# 4. Leave-one-unit-out (Table 5-12)
print("\nLeave-one-business-unit-out:")
for u in sorted(df.unit_id.unique()):
    b = fit(df[df.unit_id != u]); print(f"   excluding unit {u}: FA->DC = {b[0]:.3f}, RR->DC = {b[1]:.3f}")

# 5. Outlier treatments (Table D-4)
w = df.copy()
for c in ["zFA", "zRR", "zDC"]:
    lo, hi = w[c].quantile([0.01, 0.99]); w[c] = w[c].clip(lo, hi)
print(f"\nWinsorised 1 %: FA->DC = {fit(w)[0]:.3f}, RR->DC = {fit(w)[1]:.3f}")
nc = df[df.flag_cooks == 0]
print(f"Excluding Cook's D > 4/n (N = {len(nc)}): FA->DC = {fit(nc)[0]:.3f}, RR->DC = {fit(nc)[1]:.3f}")

# 6. Regional subgroups (Table D-3)
print("\nIndirect effect by region:")
for reg, d in df.groupby("region"):
    a = smf.ols("zRR ~ zFA", data=d).fit().params["zFA"]; b = smf.ols("zDC ~ zFA + zRR", data=d).fit().params["zRR"]
    print(f"   {reg:8s} (n={len(d):4d}): FA->RR->DC = {a*b:.3f}")

# 7. Wild-cluster bootstrap by business unit (9 clusters) for RR -> DC
X = sm.add_constant(df[["zFA", "zRR"]]); y = df.zDC
full = sm.OLS(y, X).fit(); resid = full.resid; fitted = full.fittedvalues
units = df.unit_id.values; uniq = np.unique(units)
tb = []
for _ in range(999):
    sign = dict(zip(uniq, rng.choice([-1, 1], len(uniq))))
    ystar = fitted + resid * np.array([sign[u] for u in units])
    r = sm.OLS(ystar, X).fit(cov_type="cluster", cov_kwds={"groups": units})
    tb.append((r.params["zRR"] - full.params["zRR"]) / r.bse["zRR"])
t_obs = sm.OLS(y, X).fit(cov_type="cluster", cov_kwds={"groups": units}).tvalues["zRR"]
print(f"\nWild-cluster bootstrap (9 units, Rademacher): RR->DC t = {t_obs:.2f}, bootstrap p = {(np.abs(tb) >= abs(t_obs)).mean():.3f}")
