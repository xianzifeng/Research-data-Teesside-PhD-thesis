"""
06_moderation_interaction.py  (thesis Section 4.6.1 item 5; Section 5.5.2-5.5.3; Tables 5-6, 5-7, 5-13)
  (1) H10: PR1 x PR2 interaction on FA - hierarchical regression (Model 1 main
      effects, Model 2 + interaction), simple slopes at -1 SD / mean / +1 SD of
      PR2 and the Johnson-Neyman region of significance.
  (2) H8 / H9: moderated mediation by DA and ET - group-wise a, b, c' paths with
      employee-cluster bootstrap intervals and a bootstrap difference test of the
      indirect effects, plus interaction-term tests (FA x M -> RR, RR x M -> DC).
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

from utils import OUT, CONTROLS, banner, cluster_bootstrap_indices

B = 2000
df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv"
banner("06 Moderation and interaction tests (H8-H10)", df)
cl = {"cov_type": "cluster", "cov_kwds": {"groups": df["employee_id"]}}
ctrl = " + ".join(CONTROLS)

# ---------------- (1) PR1 x PR2 -> FA (Table 5-13) ----------------------------
m1 = smf.ols(f"zFA ~ zPR1 + zPR2 + {ctrl}", data=df).fit(**cl)
m2 = smf.ols(f"zFA ~ zPR1 * zPR2 + {ctrl}", data=df).fit(**cl)
print(f"Model 1 R2 = {m1.rsquared:.3f};  Model 2 R2 = {m2.rsquared:.3f};  dR2 = {m2.rsquared - m1.rsquared:.3f}")
print(f"PR1 -> FA b = {m2.params['zPR1']:.3f};  PR2 -> FA b = {m2.params['zPR2']:.3f};  PR1 x PR2 b = {m2.params['zPR1:zPR2']:.3f} (p = {m2.pvalues['zPR1:zPR2']:.3g})")
b1, b3 = m2.params["zPR1"], m2.params["zPR1:zPR2"]
V = m2.cov_params()
print("Simple slopes of PR1 -> FA:")
for lab, zv in [("Low PR2 (-1 SD)", -1), ("Moderate (mean)", 0), ("High PR2 (+1 SD)", 1)]:
    slope = b1 + b3 * zv
    se = np.sqrt(V.loc["zPR1", "zPR1"] + zv ** 2 * V.loc["zPR1:zPR2", "zPR1:zPR2"] + 2 * zv * V.loc["zPR1", "zPR1:zPR2"])
    t = slope / se
    pr2_raw = df.PR2_perc.mean() + zv * df.PR2_perc.std()
    print(f"   {lab:18s} (PR2 perception ~ {pr2_raw:.1f}): b = {slope:.3f}, t = {t:.2f}, p = {2*stats.t.sf(abs(t), m2.df_resid):.3g}")
# Johnson-Neyman: values of zPR2 where the PR1 slope stops being significant
tcrit = stats.t.ppf(0.975, m2.df_resid)
a_ = b3 ** 2 - tcrit ** 2 * V.loc["zPR1:zPR2", "zPR1:zPR2"]
b_ = 2 * (b1 * b3 - tcrit ** 2 * V.loc["zPR1", "zPR1:zPR2"])
c_ = b1 ** 2 - tcrit ** 2 * V.loc["zPR1", "zPR1"]
disc = b_ ** 2 - 4 * a_ * c_
if disc > 0:
    roots = sorted(((-b_ - np.sqrt(disc)) / (2 * a_), (-b_ + np.sqrt(disc)) / (2 * a_)))
    raw = [df.PR2_perc.mean() + r * df.PR2_perc.std() for r in roots]
    print(f"Johnson-Neyman boundaries (zPR2): {roots[0]:.2f}, {roots[1]:.2f}  -> PR2 perception scale ~ {raw[0]:.1f}, {raw[1]:.1f}  (thesis: significant up to PR2 = 4.9)")
else:
    print("Johnson-Neyman: PR1 slope significant across the observed PR2 range")

# ---------------- (2) Moderated mediation by DA and ET (Tables 5-6, 5-7) ----------
def abc(d):
    a = sm.OLS(d.zRR, sm.add_constant(d[["zFA"]])).fit().params["zFA"]
    o = sm.OLS(d.zDC, sm.add_constant(d[["zFA", "zRR"]])).fit().params
    return np.array([a, o["zRR"], o["zFA"], a * o["zRR"]])

rng = np.random.default_rng(7)
rows = []
for mod in ["et", "da"]:
    print(f"\nModerator {mod.upper()}:")
    est, boots = {}, {}
    for g in [0, 1]:
        d = df[df[mod] == g]
        est[g] = abc(d)
        cls = d.employee_id.values
        boots[g] = np.array([abc(d.iloc[cluster_bootstrap_indices(rng, cls)]) for _ in range(B)])
        lo, hi = np.percentile(boots[g], [2.5, 97.5], axis=0)
        print(f"   {mod.upper()}={g} (N={len(d):4d}): a = {est[g][0]:.3f} [{lo[0]:.3f}, {hi[0]:.3f}]  b = {est[g][1]:.3f} [{lo[1]:.3f}, {hi[1]:.3f}]  "
              f"c' = {est[g][2]:.3f} [{lo[2]:.3f}, {hi[2]:.3f}]  indirect = {est[g][3]:.3f} [{lo[3]:.3f}, {hi[3]:.3f}]")
        rows.append({"moderator": mod, "group": g, "N": len(d), "a": est[g][0], "b": est[g][1], "c_prime": est[g][2],
                     "indirect": est[g][3], "ind_lo": lo[3], "ind_hi": hi[3]})
    diff = boots[1][:, 3] - boots[0][:, 3]
    print(f"   Difference in indirect effects ({mod.upper()}=1 - {mod.upper()}=0): {est[1][3]-est[0][3]:.3f} "
          f"[{np.percentile(diff, 2.5):.3f}, {np.percentile(diff, 97.5):.3f}]")
    ia = smf.ols(f"zRR ~ zFA * {mod} + {ctrl}", data=df).fit(**cl)
    ib = smf.ols(f"zDC ~ zFA + zRR * {mod} + {ctrl}", data=df).fit(**cl)
    print(f"   Interaction tests: FA x {mod.upper()} -> RR b = {ia.params[f'zFA:{mod}']:.3f} (p = {ia.pvalues[f'zFA:{mod}']:.3f});  "
          f"RR x {mod.upper()} -> DC b = {ib.params[f'zRR:{mod}']:.3f} (p = {ib.pvalues[f'zRR:{mod}']:.3f})")
pd.DataFrame(rows).to_csv(OUT / "06_moderated_mediation.csv", index=False)
print("\nThesis reference: H8 (DA) and H9 (ET) not supported; H10 PR1 x PR2 = -0.124 (p = 0.004), slopes 0.412 / 0.298 / 0.198")
