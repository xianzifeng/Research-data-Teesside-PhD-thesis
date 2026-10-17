"""
04_bootstrap_mediation.py  (thesis Appendix C.3; Section 4.6.1; Tables 5-4, 5-16; Section 5.3.1)
Python cross-check of the mediation estimates.  The intervals reported in
Chapter 5 are bias-corrected bootstrap intervals from Mplus
(mplus/03_mediation_bootstrap.inp, 5,000 draws).  Here we
  (a) run statsmodels' Mediation() as in Appendix C.3 (its interval is a plain
      percentile interval over independent resampling), and
  (b) run a CLUSTER bootstrap that resamples whole employees - preserving the
      within-person dependence of the six waves - and report both percentile
      and bias-corrected (BC) intervals, for the FA -> RR -> DC (H7) and
      PR -> RR -> DC (H6) pathways.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.mediation import Mediation

from utils import OUT, banner, cluster_bootstrap_indices, bc_interval

N_BOOT = 5000
df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv"
banner("04 Bootstrap mediation", df)

# (a) statsmodels Mediation, exactly as printed in Appendix C.3 ------------------
outcome_model = sm.OLS.from_formula("zDC ~ zFA + zRR", data=df)
mediator_model = sm.OLS.from_formula("zRR ~ zFA", data=df)
med = Mediation(outcome_model, mediator_model, "zFA", "zRR")
res = med.fit(n_rep=1000, method="bootstrap")       # 1,000 for speed; 5,000 in the thesis
print("statsmodels Mediation (independent resampling, percentile CI):")
print(res.summary().loc[["ACME (average)", "ADE (average)", "Total effect", "Prop. mediated (average)"]].round(3))

# (b) Cluster (employee) bootstrap for both pathways ---------------------------------
def paths(d, x):
    a = sm.OLS(d.zRR, sm.add_constant(d[[x]])).fit().params[x]
    o = sm.OLS(d.zDC, sm.add_constant(d[[x, "zRR"]])).fit().params
    return a * o["zRR"], o[x], a * o["zRR"] + o[x]            # indirect, direct, total

rng = np.random.default_rng(42)
clusters = df.employee_id.values
rows = []
for label, x in [("FA -> RR -> DC (H7)", "zFA"), ("PR -> RR -> DC (H6)", "zPR1")]:
    ind, dirc, tot = paths(df, x)
    boot = np.array([paths(df.iloc[cluster_bootstrap_indices(rng, clusters)], x)[0] for _ in range(N_BOOT)])
    pct = np.percentile(boot, [2.5, 97.5]); bc = bc_interval(boot, ind)
    print(f"\n{label}:  indirect = {ind:.3f}  direct = {dirc:.3f}  total = {tot:.3f}  PM = {ind/tot:.3f}")
    print(f"   percentile 95% CI [{pct[0]:.3f}, {pct[1]:.3f}]   bias-corrected 95% CI [{bc[0]:.3f}, {bc[1]:.3f}]   (B = {N_BOOT}, employee clusters)")
    rows.append({"pathway": label, "indirect": ind, "direct": dirc, "total": tot, "PM": ind / tot,
                 "pct_lo": pct[0], "pct_hi": pct[1], "bc_lo": bc[0], "bc_hi": bc[1], "B": N_BOOT})
pd.DataFrame(rows).to_csv(OUT / "04_mediation_bootstrap.csv", index=False)
print("\nThesis reference (Table 5-4 / 5-16): H7 indirect 0.180 [0.157, 0.205], direct 0.267, total 0.448; H6 indirect 0.071 [0.052, 0.093], total 0.205")
