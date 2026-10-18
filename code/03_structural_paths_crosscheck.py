"""
03_structural_paths_crosscheck.py  (thesis Appendix C.2; Table 5-3; Section 4.6.1, 4.9.3)
The structural model of Table 5-3 was estimated as a latent-variable SEM in
Mplus 8.8 (mplus/02_sem_structural_mlr.inp).  This script is the composite-score
cross-check referred to in Appendix D.2: path regressions on standardised
composites with cluster-robust (employee) standard errors (Arellano, 1987),
plus the business-unit-year fixed-effects specification of Table 5-11.
Composite-based coefficients differ somewhat from the latent SEM values.
"""
import pandas as pd
import statsmodels.formula.api as smf

from utils import OUT, CONTROLS, banner

df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv"
banner("03 Structural paths (composite cross-check of Table 5-3)", df)
ctrl = " + ".join(CONTROLS)
cl = {"cov_type": "cluster", "cov_kwds": {"groups": df["employee_id"]}}

models = {
    "FA ~ PR1 + PR2 + PR1xPR2  (H1/H3/H10)": f"zFA ~ zPR1 + zPR2 + zPR1xPR2 + {ctrl}",
    "RR ~ FA + PR  (H4, H6 a-path)":       f"zRR ~ zFA + zPR1 + zPR2 + {ctrl}",
    "DC ~ FA + RR + PR  (H2, H5, H1)":     f"zDC ~ zFA + zRR + zPR1 + zPR2 + {ctrl}",
}
rows = []
for label, f in models.items():
    res = smf.ols(f, data=df).fit(**cl)
    print(f"\n{label}   R2 = {res.rsquared:.3f}")
    for v in res.params.index:
        if v.startswith("z"):
            print(f"  {v:10s} b = {res.params[v]:7.3f}   t = {res.tvalues[v]:6.2f}   p = {res.pvalues[v]:.3g}")
            rows.append({"model": label, "term": v, "b": res.params[v], "se_cluster": res.bse[v], "t": res.tvalues[v], "p": res.pvalues[v]})
pd.DataFrame(rows).to_csv(OUT / "03_structural_paths.csv", index=False)

# Table 5-11: aggregated two-way (unit x year) fixed effects on the person-wave data
fe = smf.ols("zDC ~ zFA + zRR + C(unit_id) + C(wave)", data=df).fit(**cl)
print(f"\nTwo-way FE (unit, year):  FA->DC b = {fe.params['zFA']:.3f}   RR->DC b = {fe.params['zRR']:.3f}")
print("\nThesis reference (latent SEM, Table 5-3): FA->RR 0.399 | RR->DC 0.452 | FA->DC 0.267 | PR1->FA 0.412 | PR2->FA -0.198 | PR1xPR2 -0.124")
