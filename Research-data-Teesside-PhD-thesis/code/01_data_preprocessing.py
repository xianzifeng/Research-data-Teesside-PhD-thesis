"""
01_data_preprocessing.py  (thesis Appendix C.1; Sections 4.4.3, 4.4.4, 5.4)
Loads the person-wave panel, runs the quality checks reported in Section 5.4
(range checks, missingness, duplicates, VIF, |Z| > 3, Mahalanobis, Cook's D),
standardises the composites and writes the analysis file used by scripts 02-07.
"""
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm

from utils import load_panel, BLOCKS, CONTROLS, OUT, banner

df = load_panel()
banner("01 Data preprocessing and quality checks", df)
OUT.mkdir(exist_ok=True)
items = sum(BLOCKS.values(), [])

# --- 1. Range and structure checks (Section 5.4) -----------------------------
assert df[items].min().min() >= 1 and df[items].max().max() <= 5, "Likert items outside [1, 5]"
dups = df.duplicated(subset=["employee_id", "wave"]).sum()
print(f"Duplicate employee-wave rows: {dups}")
print(f"Panel balance: {df.groupby('employee_id').wave.nunique().value_counts().to_dict()} (waves per employee)")

# --- 2. Missingness (unit non-response 4-12 % per wave in the real data; FIML in Mplus) --
missing_pct = (df[items].isnull().mean() * 100).round(2)
print(f"Item non-response: max {missing_pct.max():.2f} %  (thesis: below 3 %)")
df = df.dropna(thresh=int(len(df.columns) * 0.95))        # drop rows missing > 5 % of fields

# --- 3. Standardised composites used in the Python cross-checks --------------
scaler = StandardScaler()
for c in ["FA", "RR", "DC", "PR1", "PR2"]:
    df[f"z{c}"] = scaler.fit_transform(df[[c]])
df["zPR1xPR2"] = df.zPR1 * df.zPR2

# --- 4. Multicollinearity: VIF of the structural predictors (Section 4.6.1; 5.4) ----
X = sm.add_constant(df[["zFA", "zRR", "zPR1", "zPR2"] + CONTROLS])
vif = pd.DataFrame({"Variable": X.columns,
                    "VIF": [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]})
vif = vif[vif.Variable != "const"]
print("VIF (threshold < 3 a priori; < 10 conventional):")
print(vif.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

# --- 5. Outlier diagnostics (Sections 4.4.3, 5.4) -----------------------------
z = np.abs(stats.zscore(df[["FA", "RR", "DC"]]))
print(f"Observations with |Z| > 3 on any composite: {(z > 3).any(axis=1).sum()}")
Xc = df[["FA", "RR", "DC"] + items]
cov_inv = np.linalg.pinv(np.cov(Xc.values.T))
diff = Xc.values - Xc.values.mean(axis=0)
d2 = np.einsum("ij,jk,ik->i", diff, cov_inv, diff)
crit = stats.chi2.ppf(0.999, df=Xc.shape[1])
print(f"Mahalanobis D2 > chi2(0.001, {Xc.shape[1]}) = {crit:.1f}: {(d2 > crit).sum()} cases")
ols = sm.OLS(df.zDC, sm.add_constant(df[["zFA", "zRR"]])).fit()
cooks = ols.get_influence().cooks_distance[0]
flag = cooks > 4 / len(df)
print(f"Cook's D > 4/n: {flag.sum()} observations ({flag.mean()*100:.2f} %)")
df["flag_cooks"] = flag.astype(int)

# --- 6. Save analysis file -----------------------------------------------------
df.to_csv(OUT / "analysis_panel.csv", index=False)
vif.to_csv(OUT / "01_vif.csv", index=False)
print(f"\nSaved {OUT / 'analysis_panel.csv'}  ({len(df):,} rows, {df.shape[1]} columns)")
