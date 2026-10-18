"""
00_generate_synthetic_data.py
=============================
Generates the SYNTHETIC replication dataset distributed with this repository.

The real Huafeng panel (300 employees x 6 annual waves, 2020-2025; N = 1,800
person-waves nested in 9 business units) cannot be released because of the
organisational confidentiality agreements described in Appendix H.2 of the
thesis.  This script therefore simulates a dataset that reproduces:

  * the panel STRUCTURE (ids, waves, business units, DA/ET unit attributes);
  * the 33-item measurement model of Table 4-4 (FA 5 items; RR 4 x 4 items;
    DC 4 x 3 items) plus the 6 Policy Regulation perception items (A1-A6) and
    the archival PR1/PR2 composites of Section 4.5.1;
  * the construct means / SDs / correlations of Table 5-1;
  * the standardised structural paths of Table 5-3, Table 5-13 and Section 5.7
    (PR1 -> FA 0.412, PR2 -> FA -0.198, PR1xPR2 -> FA -0.124, FA -> RR 0.399,
    RR -> DC 0.452, FA -> DC 0.267, PR -> DC 0.134, PR -> RR ~0.157);
  * within-employee clustering (ICC: FA 0.42, RR 0.38, DC 0.29; Section 4.9.3);
  * the instruments IV1/IV2 of Section 4.6.1, the four control variables of
    Section 4.5.6 and the three-item job-satisfaction marker of Section 5.4.

NOTHING in the synthetic file is a real observation.  Point estimates obtained
from it will be close to, but not identical with, the values in the thesis;
the purpose is to let every script in /code and every Mplus input in /mplus
run end-to-end.

Usage:  python code/00_generate_synthetic_data.py  [--seed 20251017]
Output: data/synthetic/synthetic_panel_person_wave.csv
        data/synthetic/synthetic_unit_year.csv
        data/synthetic/generation_parameters.json
        data/synthetic/generation_log.txt
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "synthetic"

# --------------------------------------------------------------------------
# Target parameters taken from the thesis (all cited to table / section)
# --------------------------------------------------------------------------
TARGETS = {
    "panel": {"n_employees": 300, "n_waves": 6, "first_year": 2020, "n_units": 9,
              "n_managerial": 180, "n_technical": 120},
    "composite_mean_sd": {            # Table 5-1
        "FA": (3.146, 1.183), "RR": (3.110, 0.725), "DC": (3.154, 0.814)},
    "composite_corr": {"FA_RR": 0.399, "FA_DC": 0.448, "RR_DC": 0.558},  # Table 5-1
    "paths_std": {                    # Table 5-3, Table 5-13, Section 5.7
        "PR1->FA": 0.412, "PR2->FA": -0.198, "PR1xPR2->FA": -0.124,
        "FA->RR": 0.399, "PR->RR": 0.157, "RR->DC": 0.452, "FA->DC": 0.267, "PR->DC": 0.134},
    "icc": {"FA": 0.42, "RR": 0.38, "DC": 0.29},                           # Section 4.9.3
    "loadings": {                     # Table 4-4 (standardised CFA loadings)
        "FA":  [0.786, 0.779, 0.789, 0.777, 0.786],
        "DRR": [0.777, 0.768, 0.793, 0.787],
        "TRR": [0.780, 0.778, 0.777, 0.788],
        "ORR": [0.773, 0.786, 0.763, 0.747],
        "CRR": [0.763, 0.777, 0.776, 0.765],
        "EV":  [0.809, 0.820, 0.822],
        "PSV": [0.818, 0.803, 0.815],
        "CV":  [0.816, 0.809, 0.810],
        "DIV": [0.823, 0.798, 0.805],
        "PR1p": [0.78, 0.81, 0.75],   # Appendix A items A1-A3 (EFA loadings, Table D-1)
        "PR2p": [0.72, 0.69, 0.74],   # Appendix A items A4-A6
    },
    "pr2_scale": {"mean": 3.4, "sd": 1.3},                                 # Table 5-13 (PR2 levels 2.1 / 3.4 / 4.7)
    "unit_attributes": {              # Tables 5-6 / 5-7 group sizes (person-waves)
        "ET1_share": 501 / 1800, "DA1_share": 708 / 1800},
}

# 9 business units of Table 4-2 (names are generic placeholders, not real unit names)
UNITS = pd.DataFrame({
    "unit_id": range(1, 10),
    "unit_label": ["SUB-A", "SUB-B", "SUB-C", "SUB-D", "SUB-E", "SUB-F", "SUB-G", "SUB-H", "SUB-I"],
    # ET = 1 market-facing subsidiary ; DA = 1 value-added / processed data
    "et": [1, 1, 1, 0, 0, 0, 0, 0, 0],
    "da": [1, 1, 1, 1, 0, 0, 0, 0, 0],
    "region": ["East", "East", "East", "Central", "Central", "Central", "West", "West", "East"],
    # employees per unit (sums to 300); chosen so ET=1 ~ 28 % and DA=1 ~ 39 % of person-waves
    "n_emp": [28, 28, 28, 34, 36, 36, 36, 36, 38],
})
assert UNITS.n_emp.sum() == 300


def latent_with_icc(rng, emp_index, icc, n):
    """Standard-normal latent residual with a given employee-level ICC."""
    u = rng.normal(0, np.sqrt(icc), emp_index.max() + 1)[emp_index]
    e = rng.normal(0, np.sqrt(1 - icc), n)
    return u + e


def items_from_latent(rng, eta, loadings, mean, sd_item):
    """Reflective indicators on a 1-5 Likert scale from a standardised latent score."""
    cols = []
    for lam in loadings:
        z = lam * eta + np.sqrt(1 - lam ** 2) * rng.normal(0, 1, len(eta))
        y = np.clip(np.round(mean + sd_item * z), 1, 5).astype(int)
        cols.append(y)
    return np.column_stack(cols)


def generate(seed=20251017):
    rng = np.random.default_rng(seed)
    P = TARGETS["panel"]
    n_emp, n_waves = P["n_employees"], P["n_waves"]

    # ---- employees -------------------------------------------------------
    emp = pd.DataFrame({"employee_id": [f"E{i:03d}" for i in range(1, n_emp + 1)]})
    emp["unit_id"] = np.repeat(UNITS.unit_id.values, UNITS.n_emp.values)
    role = np.array(["managerial"] * P["n_managerial"] + ["technical"] * P["n_technical"])
    rng.shuffle(role)
    emp["role"] = role
    emp["tenure_2020"] = rng.integers(1, 25, n_emp)
    emp = emp.merge(UNITS[["unit_id", "unit_label", "et", "da", "region"]], on="unit_id")

    # ---- unit-year archival variables (54 unit-years) --------------------
    years = np.arange(P["first_year"], P["first_year"] + n_waves)
    uy = UNITS[["unit_id", "unit_label", "et", "da", "region"]].merge(
        pd.DataFrame({"year": years}), how="cross")
    uy["wave"] = uy["year"] - P["first_year"] + 1
    unit_eff = rng.normal(0, 0.6, 9)
    trend = (uy["wave"] - 3.5) * 0.18                       # policy support rises over 2020-2025
    # PR1 (policy support intensity) and PR2 (compliance constraints) as archival composites
    uy["pr1_pilots"] = np.clip(np.round(2 + 0.6 * trend + unit_eff[uy.unit_id - 1] + rng.normal(0, 0.8, len(uy))), 0, 8).astype(int)
    uy["pr1_log_subsidy"] = np.log1p(np.exp(5.5 + 0.4 * trend + 0.5 * unit_eff[uy.unit_id - 1] + rng.normal(0, 0.5, len(uy))))
    uy["pr1_tax_incentive"] = np.clip(np.round(3 + 0.3 * trend + 0.5 * unit_eff[uy.unit_id - 1] + rng.normal(0, 0.9, len(uy))), 1, 5).astype(int)
    uy["pr2_n_requirements"] = np.clip(np.round(6 + 0.5 * trend + rng.normal(0, 1.5, len(uy))), 2, 12).astype(int)
    uy["pr2_penalty_severity"] = np.clip(np.round(3.2 + 0.2 * trend + rng.normal(0, 1.0, len(uy))), 1, 5).astype(int)
    uy["pr2_free_data_scope"] = np.clip(np.round(3.4 + rng.normal(0, 1.0, len(uy))), 1, 5).astype(int)
    # Controls (Section 4.5.6)
    uy["ctrl_log_employees"] = np.log1p(np.exp(5.2 + 0.4 * unit_eff[uy.unit_id - 1] + 0.05 * (uy.wave - 1)))
    uy["ctrl_firm_age"] = (10 + uy.unit_id * 1.5 + (uy.year - 2020)).astype(int)
    uy["ctrl_n_competitors"] = np.clip(np.round(6 + rng.normal(0, 2.5, len(uy))), 1, 15).astype(int)
    uy["ctrl_log_gdp_pc"] = np.log1p(np.exp(11.2 + uy.region.map({"East": 0.35, "Central": 0.0, "West": -0.3}) + rng.normal(0, 0.05, len(uy))))
    # Instruments (Section 4.6.1): relevant to PR1 by construction
    uy["iv1_adj_province_support"] = 0.7 * (uy.pr1_pilots - uy.pr1_pilots.mean()) / uy.pr1_pilots.std() + rng.normal(0, 0.7, len(uy))
    uy["iv2_obs_station_density"] = 0.5 * unit_eff[uy.unit_id - 1] + rng.normal(0, 0.8, len(uy)) + 4.0
    # Archival outcome indicators used as 4th / 5th composite indicators (Section 4.5.4 / assumption 2 of the revision notes)
    uy["arch_dc_revenue_growth_pct"] = np.round(8 + 3 * trend + 2 * unit_eff[uy.unit_id - 1] + rng.normal(0, 4, len(uy)), 1)
    uy["arch_dc_direct_sales_share_pct"] = np.round(np.clip(14.5 + 1.5 * trend + 3 * uy.da + rng.normal(0, 3, len(uy)), 2, 40), 1)

    # Standardised PR composites at unit-year level
    z = lambda s: (s - s.mean()) / s.std()
    uy["PR1"] = z(z(uy.pr1_pilots) + z(uy.pr1_log_subsidy) + z(uy.pr1_tax_incentive))
    uy["PR2"] = z(z(uy.pr2_n_requirements) + z(uy.pr2_penalty_severity) + z(uy.pr2_free_data_scope))
    uy["PR"] = z(uy.PR1 - uy.PR2)   # net policy regulation score used for PR -> RR / PR -> DC paths

    # ---- person-wave panel -----------------------------------------------
    pw = emp.merge(pd.DataFrame({"wave": np.arange(1, n_waves + 1)}), how="cross")
    pw["year"] = pw.wave + P["first_year"] - 1
    pw = pw.merge(uy.drop(columns=["unit_label", "et", "da", "region"]), on=["unit_id", "wave", "year"])
    pw = pw.sort_values(["employee_id", "wave"]).reset_index(drop=True)
    n = len(pw)
    emp_idx = pd.factorize(pw.employee_id)[0]

    T, I = TARGETS["paths_std"], TARGETS["icc"]
    pr1, pr2 = pw.PR1.values, pw.PR2.values
    # Latent structural model (standardised, recursive)
    fa_lat = T["PR1->FA"] * pr1 + T["PR2->FA"] * pr2 + T["PR1xPR2->FA"] * (pr1 * pr2)
    fa = fa_lat + np.sqrt(max(1 - np.var(fa_lat), 0.3)) * latent_with_icc(rng, emp_idx, min(I["FA"] * 1.25, 0.9), n)
    fa = (fa - fa.mean()) / fa.std()
    rr_lat = (T["FA->RR"] - 0.04) * fa + T["PR->RR"] * pw.PR.values  # -0.04 offsets covariance added via PR
    rr = rr_lat + np.sqrt(max(1 - np.var(rr_lat), 0.3)) * latent_with_icc(rng, emp_idx, min(I["RR"] * 1.25, 0.9), n)
    rr = (rr - rr.mean()) / rr.std()
    dc_lat = T["RR->DC"] * rr + T["FA->DC"] * fa + T["PR->DC"] * pw.PR.values
    dc = dc_lat + np.sqrt(max(1 - np.var(dc_lat), 0.3)) * latent_with_icc(rng, emp_idx, min(I["DC"] * 0.95, 0.9), n)
    dc = (dc - dc.mean()) / dc.std()
    # RR and DC sub-dimensions: second-order structure (first-order factors load ~0.85 on the construct)
    sub = lambda eta: (0.85 * eta + np.sqrt(1 - 0.85 ** 2) * rng.normal(0, 1, n))

    L = TARGETS["loadings"]
    M = TARGETS["composite_mean_sd"]
    blocks = {}
    # item SDs chosen so that the composite SDs approximate Table 5-1 after rounding
    blocks["fa"] = items_from_latent(rng, fa, L["FA"], M["FA"][0], 1.62)
    for d in ["DRR", "TRR", "ORR", "CRR"]:
        blocks[d.lower()] = items_from_latent(rng, sub(rr), L[d], M["RR"][0], 1.10)
    for d in ["EV", "PSV", "CV", "DIV"]:
        blocks[d.lower()] = items_from_latent(rng, sub(dc), L[d], M["DC"][0], 1.22)
    blocks["a_pr1"] = items_from_latent(rng, pr1, L["PR1p"], 3.2, 1.1)
    blocks["a_pr2"] = items_from_latent(rng, pr2, L["PR2p"], TARGETS["pr2_scale"]["mean"], TARGETS["pr2_scale"]["sd"])
    # Marker variable: job satisfaction (3 MSQ-type items), unrelated to the model
    js = rng.normal(0, 1, n)
    blocks["js"] = items_from_latent(rng, js, [0.80, 0.78, 0.82], 3.5, 1.0)

    names = {"fa": [f"fa{i}" for i in range(1, 6)],
             "drr": [f"drr{i}" for i in range(1, 5)], "trr": [f"trr{i}" for i in range(1, 5)],
             "orr": [f"orr{i}" for i in range(1, 5)], "crr": [f"crr{i}" for i in range(1, 5)],
             "ev": [f"ev{i}" for i in range(1, 4)], "psv": [f"psv{i}" for i in range(1, 4)],
             "cv": [f"cv{i}" for i in range(1, 4)], "div": [f"div{i}" for i in range(1, 4)],
             "a_pr1": ["a1", "a2", "a3"], "a_pr2": ["a4", "a5", "a6"], "js": ["js1", "js2", "js3"]}
    for k, cols in names.items():
        pw[cols] = blocks[k]

    # Composite scores (mean of items), as used by the Python cross-checks
    pw["FA"] = pw[names["fa"]].mean(axis=1)
    for d in ["drr", "trr", "orr", "crr", "ev", "psv", "cv", "div"]:
        pw[d.upper()] = pw[names[d]].mean(axis=1)
    pw["RR"] = pw[["DRR", "TRR", "ORR", "CRR"]].mean(axis=1)
    pw["DC"] = pw[["EV", "PSV", "CV", "DIV"]].mean(axis=1)
    pw["PR1_perc"] = pw[["a1", "a2", "a3"]].mean(axis=1)
    pw["PR2_perc"] = pw[["a4", "a5", "a6"]].mean(axis=1)
    pw["JS_marker"] = pw[["js1", "js2", "js3"]].mean(axis=1)
    pw["PR1xPR2"] = pw.PR1 * pw.PR2

    # Order columns
    front = ["employee_id", "unit_id", "unit_label", "wave", "year", "role", "tenure_2020", "et", "da", "region"]
    items = sum(names.values(), [])
    comps = ["FA", "DRR", "TRR", "ORR", "CRR", "RR", "EV", "PSV", "CV", "DIV", "DC", "PR1_perc", "PR2_perc", "JS_marker"]
    arch = ["PR1", "PR2", "PR", "PR1xPR2"] + [c for c in uy.columns if c.startswith(("pr1_", "pr2_", "ctrl_", "iv", "arch_"))]
    pw = pw[front + items + comps + arch]
    return pw, uy, emp


def report(pw):
    lines = []
    M = TARGETS["composite_mean_sd"]
    for c in ["FA", "RR", "DC"]:
        lines.append(f"{c}: mean {pw[c].mean():.3f} (target {M[c][0]}), sd {pw[c].std():.3f} (target {M[c][1]})")
    corr = pw[["FA", "RR", "DC"]].corr()
    lines.append(f"corr FA-RR {corr.loc['FA','RR']:.3f} (0.399), FA-DC {corr.loc['FA','DC']:.3f} (0.448), RR-DC {corr.loc['RR','DC']:.3f} (0.558)")

    def alpha(df):
        k = df.shape[1]; v = df.var(ddof=1)
        return k / (k - 1) * (1 - v.sum() / df.sum(axis=1).var(ddof=1))
    for blk, cols in {"FA": ["fa1", "fa2", "fa3", "fa4", "fa5"], "DRR": ["drr1", "drr2", "drr3", "drr4"],
                      "EV": ["ev1", "ev2", "ev3"], "PR2 (A4-A6)": ["a4", "a5", "a6"]}.items():
        lines.append(f"alpha {blk}: {alpha(pw[cols]):.3f}")

    def icc1(y, g):
        d = pd.DataFrame({"y": y, "g": g}); m = d.groupby("g").y.agg(["mean", "count"])
        msb = (m["count"] * (m["mean"] - d.y.mean()) ** 2).sum() / (len(m) - 1)
        msw = ((d.y - d.g.map(m["mean"])) ** 2).sum() / (len(d) - len(m))
        k = m["count"].mean(); return (msb - msw) / (msb + (k - 1) * msw)
    for c in ["FA", "RR", "DC"]:
        lines.append(f"ICC(1) {c}: {icc1(pw[c], pw.employee_id):.3f} (target {TARGETS['icc'][c]})")
    lines.append(f"ET=1 person-waves: {int(pw.et.sum())} (thesis 501); DA=1: {int(pw.da.sum())} (thesis 708)")
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20251017)
    a = ap.parse_args()
    pw, uy, emp = generate(a.seed)
    OUT.mkdir(parents=True, exist_ok=True)
    pw.to_csv(OUT / "synthetic_panel_person_wave.csv", index=False)
    uy.to_csv(OUT / "synthetic_unit_year.csv", index=False)
    (OUT / "generation_parameters.json").write_text(json.dumps({"seed": a.seed, "targets": TARGETS}, indent=2))
    log = report(pw)
    (OUT / "generation_log.txt").write_text("Synthetic data generation log (seed %d)\n%s\n" % (a.seed, log))
    print(log)
