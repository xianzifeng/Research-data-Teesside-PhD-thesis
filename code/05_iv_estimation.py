"""
05_iv_estimation.py  (thesis Appendix C.4; Section 4.6.1 item 4; Table 5-3 Panel B; Section 6.5.2)
2SLS instrumental-variable regressions with employee-clustered standard errors.
Equation 1 instruments FA, equation 2 instruments RR, both with
  IV1 = adjacent-province policy support intensity (iv1_adj_province_support)
  IV2 = regional observation-station density      (iv2_obs_station_density)
and the four controls of Section 4.5.6.  First-stage F and the Sargan/Hansen J
statistic indicate relevance and internal consistency; neither establishes the
exclusion restriction, so the IV estimates are reported as assumption-dependent
robustness checks, not as causal identification (Section 6.5.2).
"""
import pandas as pd
from linearmodels.iv import IV2SLS

from utils import OUT, CONTROLS, banner

df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv"
banner("05 Instrumental-variable estimation (Table 5-3, Panel B)", df)
ctrl = " + ".join(CONTROLS)
eqs = {
    "Eq.1  DC ~ [FA ~ IV1 + IV2] + RR": f"zDC ~ 1 + [zFA ~ iv1_adj_province_support + iv2_obs_station_density] + zRR + {ctrl}",
    "Eq.2  DC ~ [RR ~ IV1 + IV2] + FA": f"zDC ~ 1 + [zRR ~ iv1_adj_province_support + iv2_obs_station_density] + zFA + {ctrl}",
}
rows = []
for label, f in eqs.items():
    res = IV2SLS.from_formula(f, data=df).fit(cov_type="clustered", clusters=df["employee_id"])
    endog = "zFA" if "[zFA" in f else "zRR"
    fs = res.first_stage.diagnostics
    print(f"\n{label}")
    print(f"   beta_IV({endog} -> DC) = {res.params[endog]:.3f}  (SE {res.std_errors[endog]:.3f}, p = {res.pvalues[endog]:.3g})")
    print(f"   first-stage F = {fs.loc[endog, 'f.stat']:.1f}  partial R2 = {fs.loc[endog, 'rsquared']:.3f}")
    print(f"   Sargan J = {res.sargan.stat:.2f} (p = {res.sargan.pval:.3f})   Wooldridge endogeneity p = {res.wooldridge_regression.pval:.3f}")
    rows.append({"equation": label, "endog": endog, "beta_iv": res.params[endog], "se": res.std_errors[endog],
                 "first_stage_F": fs.loc[endog, "f.stat"], "sargan_J": res.sargan.stat, "sargan_p": res.sargan.pval})
pd.DataFrame(rows).to_csv(OUT / "05_iv_estimates.csv", index=False)
print("\nThesis reference (composite-based, Section 6.1.3): FA->DC beta_IV 0.312 vs OLS 0.287; RR->DC 0.356 vs 0.341")
