"""
02_measurement_checks.py  (thesis Sections 4.4.1 Stage 4, 4.5.7, 4.7.1, 4.9.3, 5.4; Tables 4-3, 4-4, 5-2, D-1)
Python cross-checks of the measurement model.  The CFA loadings, CR/AVE and fit
indices reported in Tables 4-4 and 5-2 come from Mplus 8.8 (mplus/01_cfa_measurement.inp);
this script reproduces the reliability/validity statistics that do not require
a latent-variable estimator and the Wave-1 EFA / KMO / Bartlett / ICC / Harman diagnostics.
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from utils import (BLOCKS, SECOND_ORDER, OUT, banner, cronbach_alpha, standardized_loadings,
                   cr_ave, htmt, icc1, kmo_bartlett)

df = pd.read_csv(OUT / "analysis_panel.csv")
df.attrs["source"] = "analysis_panel.csv (see 01)"
banner("02 Measurement model checks", df)

# --- 1. Reliability and convergent validity (Table 5-2) ----------------------
rows = []
for blk, cols in BLOCKS.items():
    lam = standardized_loadings(df[cols])
    cr, ave = cr_ave(lam)
    rows.append({"Construct": blk, "Items": len(cols), "Cronbach_alpha": cronbach_alpha(df[cols]),
                 "CR": cr, "AVE": ave, "Loading_min": lam.min(), "Loading_max": lam.max()})
rel = pd.DataFrame(rows)
print(rel.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
rel.to_csv(OUT / "02_reliability_validity.csv", index=False)

# --- 2. Discriminant validity: Fornell-Larcker and HTMT (Section 4.5.7) -------
first_order = [b for b in BLOCKS if not b.startswith("PR")]
scores = pd.DataFrame({b: df[BLOCKS[b]].mean(axis=1) for b in first_order})
R = scores.corr()
sqrt_ave = {r["Construct"]: np.sqrt(r["AVE"]) for _, r in rel.iterrows()}
viol = [(a, b) for a in first_order for b in first_order if a < b and R.loc[a, b] > min(sqrt_ave[a], sqrt_ave[b])]
print(f"\nFornell-Larcker violations (r > sqrt(AVE)): {len(viol)}  {viol if viol else ''}")
H = pd.DataFrame(index=first_order, columns=first_order, dtype=float)
for a in first_order:
    for b in first_order:
        if a < b:
            H.loc[a, b] = htmt(df, BLOCKS[a], BLOCKS[b])
print(f"HTMT range: {np.nanmin(H.values):.3f} - {np.nanmax(H.values):.3f}  (threshold 0.85, Henseler et al., 2015)")
H.to_csv(OUT / "02_htmt_matrix.csv")

# --- 3. Wave-1 EFA (Table D-1): eigenvalues, KMO, Bartlett ---------------------
w1 = df[df.wave == 1]
efa_items = sum([BLOCKS[b] for b in ["PR1_perc", "PR2_perc", "FA", "DRR", "TRR", "ORR", "CRR", "EV", "PSV", "CV", "DIV"]], [])
kmo, chi2, dof, p = kmo_bartlett(w1[efa_items])
eig = PCA().fit((w1[efa_items] - w1[efa_items].mean()) / w1[efa_items].std()).explained_variance_
print(f"\nWave-1 EFA (n = {len(w1)}, {len(efa_items)} items): KMO = {kmo:.3f}; Bartlett chi2({int(dof)}) = {chi2:.1f}, p = {p:.3g}")
print("Eigenvalues > 1:", np.round(eig[eig > 1], 2), f"; cumulative variance of first 4: {eig[:4].sum()/eig.sum()*100:.1f} %")
# Parallel analysis (recommended in the revision notes, comment C769)
rng = np.random.default_rng(1)
rand_eig = np.mean([PCA().fit(rng.normal(size=w1[efa_items].shape)).explained_variance_ for _ in range(100)], axis=0)
print(f"Parallel analysis retains {int((eig > rand_eig).sum())} components")

# --- 4. Within-employee clustering (Section 4.9.3) -----------------------------
print()
for c in ["FA", "RR", "DC"]:
    rho = icc1(df[c], df.employee_id)
    n_eff = len(df) / (1 + (6 - 1) * rho)
    print(f"ICC(1) {c} = {rho:.3f}  ->  N_effective = {n_eff:.0f}")

# --- 5. Common method bias: Harman single factor and marker variable (Section 5.4) ---
all_items = sum([BLOCKS[b] for b in first_order], [])
pca = PCA(n_components=1).fit((df[all_items] - df[all_items].mean()) / df[all_items].std())
print(f"\nHarman single-factor: first unrotated factor explains {pca.explained_variance_ratio_[0]*100:.1f} % of variance (thesis 28.4 %)")
m = df["JS_marker"]
def partial_r(x, y, z):
    rxy, rxz, ryz = np.corrcoef(x, y)[0, 1], np.corrcoef(x, z)[0, 1], np.corrcoef(y, z)[0, 1]
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
for a, b in [("FA", "RR"), ("RR", "DC")]:
    print(f"Marker-variable check {a}-{b}: r = {np.corrcoef(df[a], df[b])[0,1]:.3f}, r_partial|JS = {partial_r(df[a], df[b], m):.3f}")
