"""
10_sem_latent_semopy.py  -- open-source latent-variable SEM equivalent of the five Mplus inputs in /mplus
=====================================================================================================
Fits the SAME models as mplus/01-05 with semopy (Python), so that the latent-variable
results of Tables 4-4, 5-2, 5-3, 5-4, 5-6, 5-7, 5-13 and 5-16 can be produced without
an Mplus licence.  Standard errors and confidence intervals come from an EMPLOYEE-CLUSTER
BOOTSTRAP (whole employees resampled with replacement), which respects the person-wave
nesting exactly as TYPE = COMPLEX / CLUSTER = emp does in Mplus.

Models
  M1  second-order CFA (PR1P, PR2P, FA, DRR/TRR/ORR/CRR -> RR, EV/PSV/CV/DIV -> DC)         [mplus/01]
  M2  structural model  FA ON PR1 PR2 ctrl ; RR ON FA PR1 PR2 ctrl ; DC ON FA RR PR1 PR2 ctrl  [mplus/02]
  M3  indirect effects  FA->RR->DC (H7), PR1->RR->DC (H6) with percentile and BC bootstrap CIs [mplus/03]
  M4  interaction       FA ON PR1 PR2 PR1xPR2 ; DC ON ... PR1xPR2 ; simple slopes             [mplus/04]
  M5  multi-group       a, b, c' and indirect by ET and by DA, bootstrap difference            [mplus/05]

Usage
  python code/10_sem_latent_semopy.py                 # synthetic data, B = 200 (quick check, ~10 min)
  python code/10_sem_latent_semopy.py --data /path/to/real_panel.csv --boot 5000 --out outputs/real_results
  python code/10_sem_latent_semopy.py --models M1 M2   # subset

The real panel never leaves the author's computer: pass its path with --data; outputs are
summary tables only (no row-level data are written).
Requires: semopy >= 2.3 (pip install semopy), pandas, numpy, scipy, statsmodels.
"""
import argparse, json, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
import semopy

from utils import BLOCKS, CONTROLS, cluster_bootstrap_indices, bc_interval, ROOT, DATA_SYNTH

# ----------------------------------------------------------------------------------------------
# Model syntax (lavaan-style, as used by semopy and directly comparable to the Mplus BY/ON lines)
# ----------------------------------------------------------------------------------------------
def meas(include_pr=False):
    lines = [f"{k} =~ " + " + ".join(v) for k, v in BLOCKS.items() if include_pr or not k.startswith("PR")]
    lines += ["RR =~ DRR + TRR + ORR + CRR", "DC =~ EV + PSV + CV + DIV"]
    return "\n".join(lines).replace("PR1_perc =~", "PR1P =~").replace("PR2_perc =~", "PR2P =~")

CTRL = " + ".join(CONTROLS)
M1 = meas(include_pr=True)
M2 = meas() + f"""
FA ~ PR1 + PR2 + {CTRL}
RR ~ FA + PR1 + PR2 + {CTRL}
DC ~ FA + RR + PR1 + PR2 + {CTRL}
"""
M4 = meas() + f"""
FA ~ PR1 + PR2 + PR1xPR2 + {CTRL}
RR ~ FA + PR1 + PR2 + {CTRL}
DC ~ FA + RR + PR1 + PR2 + PR1xPR2 + {CTRL}
"""
M5 = meas() + """
RR ~ FA
DC ~ FA + RR
"""


def fit(model_desc, data, obj="MLW"):
    m = semopy.Model(model_desc)
    m.fit(data, obj=obj)
    return m


def std_paths(m):
    """Standardised regression coefficients (STDYX) as a dict {'lhs~rhs': value}."""
    # semopy writes loadings as "indicator ~ latent"; structural paths have FA / RR / DC on the left
    ins = m.inspect(std_est=True)
    reg = ins[(ins.op == "~") & ins.lval.isin(["FA", "RR", "DC"])]
    return {f"{r.lval}~{r.rval}": float(r["Est. Std"]) for _, r in reg.iterrows()}


def std_loadings(m):
    ins = m.inspect(std_est=True)
    lo = ins[(ins.op == "~") & ~ins.lval.isin(["FA", "RR", "DC"]) | ((ins.op == "~") & ins.lval.isin(["DRR", "TRR", "ORR", "CRR", "EV", "PSV", "CV", "DIV"]))]
    return {f"{r.rval}=~{r.lval}": float(r["Est. Std"]) for _, r in lo.iterrows()}


def fit_indices(m):
    s = semopy.calc_stats(m).T.iloc[:, 0]
    return {k: float(s[k]) for k in ["chi2", "DoF", "chi2 p-value", "CFI", "TLI", "RMSEA", "AIC", "BIC"] if k in s}


def srmr(m, data):
    """Standardised root mean square residual (semopy does not report it)."""
    S = data[m.vars["observed"]].corr().values
    Sigma = m.calc_sigma()[0]
    d = np.sqrt(np.diag(Sigma)); Sig = Sigma / np.outer(d, d)
    iu = np.triu_indices_from(S)
    return float(np.sqrt(np.mean((S[iu] - Sig[iu]) ** 2)))


def prepare(df):
    d = df.copy()
    for c in ["PR1", "PR2"]:
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    d["PR1xPR2"] = d.PR1 * d.PR2
    for c in CONTROLS:
        d[c] = (d[c] - d[c].mean()) / d[c].std()
    return d


def cluster_boot(model_desc, data, extract, B, seed, obj="MLW", label=""):
    """Refit the model on B employee-cluster resamples; return array of extracted statistics."""
    rng = np.random.default_rng(seed)
    cl = data.employee_id.values
    out, fails, t0 = [], 0, time.time()
    for b in range(B):
        try:
            d = data.iloc[cluster_bootstrap_indices(rng, cl)].reset_index(drop=True)
            out.append(extract(fit(model_desc, d, obj)))
        except Exception:
            fails += 1
        if (b + 1) % 50 == 0:
            print(f"      {label} bootstrap {b+1}/{B}  ({time.time()-t0:.0f}s, {fails} failed)")
    return np.array(out), fails


def summarise(est: dict, boot: np.ndarray, keys):
    rows = []
    for i, k in enumerate(keys):
        b = boot[:, i]
        lo, hi = np.percentile(b, [2.5, 97.5]); bl, bh = bc_interval(b, est[k])
        se = b.std(ddof=1)
        rows.append({"parameter": k, "estimate": est[k], "boot_se": se, "z": est[k] / se if se > 0 else np.nan,
                     "p_boot": 2 * (1 - stats.norm.cdf(abs(est[k] / se))) if se > 0 else np.nan,
                     "pct_lo": lo, "pct_hi": hi, "bc_lo": bl, "bc_hi": bh})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=None, help="CSV with the codebook columns; default = synthetic panel")
    ap.add_argument("--boot", type=int, default=200, help="bootstrap draws (thesis: 5000)")
    ap.add_argument("--models", nargs="*", default=["M1", "M2", "M3", "M4", "M5"])
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=2026)
    a = ap.parse_args()

    path = Path(a.data) if a.data else DATA_SYNTH
    src = "REAL" if a.data else "SYNTHETIC"
    out = Path(a.out) if a.out else ROOT / "outputs" / ("sem_real" if a.data else "sem_synthetic")
    out.mkdir(parents=True, exist_ok=True)
    df = prepare(pd.read_csv(path))
    print(f"semopy {semopy.__version__} | data: {src} ({len(df)} rows, {df.employee_id.nunique()} employees) | B = {a.boot}")
    log = {"data": src, "n_rows": len(df), "n_clusters": int(df.employee_id.nunique()), "B": a.boot, "semopy": semopy.__version__}

    # ---------------- M1: CFA ----------------
    if "M1" in a.models:
        print("\n[M1] Second-order CFA (mplus/01)")
        m = fit(M1, df)
        fi = fit_indices(m); fi["SRMR"] = srmr(m, df)
        print("   fit:", {k: round(v, 3) for k, v in fi.items()})
        lo = std_loadings(m)
        rows = []
        for blk, items in BLOCKS.items():
            name = {"PR1_perc": "PR1P", "PR2_perc": "PR2P"}.get(blk, blk)
            lam = np.array([abs(lo[f"{name}=~{i}"]) for i in items])
            err = 1 - lam ** 2
            rows.append({"construct": blk, "items": len(items), "loading_min": lam.min(), "loading_max": lam.max(),
                         "CR": lam.sum() ** 2 / (lam.sum() ** 2 + err.sum()), "AVE": (lam ** 2).mean()})
        for so, dims in {"RR": ["DRR", "TRR", "ORR", "CRR"], "DC": ["EV", "PSV", "CV", "DIV"]}.items():
            lam = np.array([abs(lo[f"{so}=~{d}"]) for d in dims])
            rows.append({"construct": f"{so} (second-order)", "items": len(dims), "loading_min": lam.min(), "loading_max": lam.max(),
                         "CR": lam.sum() ** 2 / (lam.sum() ** 2 + (1 - lam ** 2).sum()), "AVE": (lam ** 2).mean()})
        tab = pd.DataFrame(rows); print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        tab.to_csv(out / "M1_cfa_reliability.csv", index=False)
        pd.Series(lo).to_csv(out / "M1_cfa_std_loadings.csv", header=["std_loading"])
        log["M1_fit"] = fi

    # ---------------- M2: structural model ----------------
    keys2 = ["FA~PR1", "FA~PR2", "RR~FA", "RR~PR1", "RR~PR2", "DC~FA", "DC~RR", "DC~PR1", "DC~PR2"]
    if "M2" in a.models or "M3" in a.models:
        print("\n[M2] Structural model (mplus/02) + [M3] indirect effects (mplus/03)")
        m = fit(M2, df)
        fi = fit_indices(m); fi["SRMR"] = srmr(m, df); print("   fit:", {k: round(v, 3) for k, v in fi.items()})
        p = std_paths(m)
        est = {k: p[k] for k in keys2}
        est["IND_FA_RR_DC"] = p["RR~FA"] * p["DC~RR"]
        est["TOT_FA_DC"] = est["IND_FA_RR_DC"] + p["DC~FA"]
        est["IND_PR1_RR_DC"] = p["RR~PR1"] * p["DC~RR"]
        est["IND_PR1_FA_DC"] = p["FA~PR1"] * p["DC~FA"]
        est["IND_PR1_FA_RR_DC"] = p["FA~PR1"] * p["RR~FA"] * p["DC~RR"]
        est["TOT_PR1_DC"] = p["DC~PR1"] + est["IND_PR1_RR_DC"] + est["IND_PR1_FA_DC"] + est["IND_PR1_FA_RR_DC"]
        keys = list(est.keys())
        def ex(mm):
            q = std_paths(mm)
            v = [q[k] for k in keys2]
            ind = q["RR~FA"] * q["DC~RR"]
            v += [ind, ind + q["DC~FA"], q["RR~PR1"] * q["DC~RR"], q["FA~PR1"] * q["DC~FA"],
                  q["FA~PR1"] * q["RR~FA"] * q["DC~RR"]]
            v += [q["DC~PR1"] + v[-3] + v[-2] + v[-1]]
            return v
        boot, fails = cluster_boot(M2, df, ex, a.boot, a.seed, label="M2/M3")
        tab = summarise(est, boot, keys)
        print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        tab.to_csv(out / "M2_M3_structural_and_indirect.csv", index=False)
        log["M2_fit"] = fi; log["M2_boot_failed"] = fails
        print("   thesis reference: PR1->FA 0.412, PR2->FA -0.198, FA->RR 0.399, RR->DC 0.452, FA->DC 0.267, PR->DC 0.134;")
        print("                     H7 indirect 0.180 [0.157, 0.205]; H6 indirect 0.071 [0.052, 0.093]")

    # ---------------- M4: PR1 x PR2 interaction ----------------
    if "M4" in a.models:
        print("\n[M4] PR1 x PR2 interaction (mplus/04)")
        m = fit(M4, df); p = std_paths(m)
        keys4 = ["FA~PR1", "FA~PR2", "FA~PR1xPR2", "DC~PR1xPR2"]
        est = {k: p[k] for k in keys4}
        est["slope_low"] = p["FA~PR1"] - p["FA~PR1xPR2"]
        est["slope_mean"] = p["FA~PR1"]
        est["slope_high"] = p["FA~PR1"] + p["FA~PR1xPR2"]
        keys = list(est.keys())
        def ex4(mm):
            q = std_paths(mm); b1, b3 = q["FA~PR1"], q["FA~PR1xPR2"]
            return [q[k] for k in keys4] + [b1 - b3, b1, b1 + b3]
        boot, fails = cluster_boot(M4, df, ex4, a.boot, a.seed + 1, label="M4")
        tab = summarise(est, boot, keys)
        print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        tab.to_csv(out / "M4_interaction.csv", index=False); log["M4_boot_failed"] = fails
        print("   thesis reference: PR1xPR2 -> FA = -0.124 (p = 0.004); slopes 0.412 / 0.298 / 0.198")

    # ---------------- M5: multi-group by ET and DA ----------------
    if "M5" in a.models:
        for g in ["et", "da"]:
            print(f"\n[M5] Multi-group by {g.upper()} (mplus/05)")
            res, boots = {}, {}
            for lvl in [0, 1]:
                d = df[df[g] == lvl].reset_index(drop=True)
                m = fit(M5, d); p = std_paths(m)
                res[lvl] = {"a": p["RR~FA"], "b": p["DC~RR"], "c_prime": p["DC~FA"], "indirect": p["RR~FA"] * p["DC~RR"], "N": len(d)}
                def ex5(mm):
                    q = std_paths(mm); return [q["RR~FA"], q["DC~RR"], q["DC~FA"], q["RR~FA"] * q["DC~RR"]]
                boots[lvl], _ = cluster_boot(M5, d, ex5, a.boot, a.seed + 2 + lvl, label=f"M5 {g}={lvl}")
            rows = []
            for lvl in [0, 1]:
                lo, hi = np.percentile(boots[lvl], [2.5, 97.5], axis=0)
                r = res[lvl]
                rows.append({"group": f"{g}={lvl}", "N": r["N"], "a": r["a"], "a_lo": lo[0], "a_hi": hi[0], "b": r["b"], "b_lo": lo[1], "b_hi": hi[1],
                             "c_prime": r["c_prime"], "c_lo": lo[2], "c_hi": hi[2], "indirect": r["indirect"], "ind_lo": lo[3], "ind_hi": hi[3]})
            nb = min(len(boots[0]), len(boots[1]))
            diff = boots[1][:nb, 3] - boots[0][:nb, 3]
            rows.append({"group": f"diff ({g}=1 - {g}=0)", "indirect": res[1]["indirect"] - res[0]["indirect"],
                         "ind_lo": np.percentile(diff, 2.5), "ind_hi": np.percentile(diff, 97.5)})
            tab = pd.DataFrame(rows); print(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
            tab.to_csv(out / f"M5_multigroup_{g}.csv", index=False)
        print("   thesis reference: ET difference -0.026 [-0.077, 0.028]; DA difference 0.016 [-0.030, 0.064]")

    (out / "run_info.json").write_text(json.dumps(log, indent=2, default=float))
    print(f"\nTables written to {out}")


if __name__ == "__main__":
    main()
