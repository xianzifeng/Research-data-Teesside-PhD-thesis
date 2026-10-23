"""
11_prepare_real_panel.py  -- convert the author's raw survey/archival files into the codebook layout
====================================================================================================
Run this ONLY on the author's own computer.  It reads the raw Excel export named in Appendix C.1
(huafeng_survey_2020_2025.xlsx) plus the unit-year archival sheet, renames columns to the names
in docs/CODEBOOK.md and writes data/real/huafeng_panel_person_wave.csv, which is ignored by git
(.gitignore) and therefore can never be pushed by accident.

EDIT THE TWO DICTIONARIES BELOW so that the left-hand side is the column name in your raw file
and the right-hand side is the codebook name.  Leave out anything you do not have; the SEM
script (10_sem_latent_semopy.py) needs at minimum:
   employee_id, wave, unit_id, et, da,
   a1-a6, fa1-fa5(*), drr1-4, trr1-4, orr1-4, crr1-4, ev1-3, psv1-3, cv1-3, div1-3,
   PR1, PR2, ctrl_log_employees, ctrl_firm_age, ctrl_n_competitors, ctrl_log_gdp_pc
(*) fa5 and drr4/trr4/orr4/crr4 are the archival / expert-scored composite indicators of the
    mixed measurement model (Table 4-4); include them as already-computed columns.

Usage:  python code/11_prepare_real_panel.py --survey ~/path/huafeng_survey_2020_2025.xlsx \
            --sheet raw_data --archival ~/path/unit_year_archival.xlsx
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "real" / "huafeng_panel_person_wave.csv"

# --- 1. survey file: raw column name -> codebook name -------------------------------------------
SURVEY_COLUMNS = {
    # identifiers
    "respondent_id": "employee_id",      # pseudonymous id, same across waves
    "survey_year": "year",               # 2020 ... 2025  (wave is derived)
    "business_unit": "unit_id",          # 1 ... 9
    "role": "role",
    # Policy Regulation perception items (Appendix A, A1-A6)
    "A1": "a1", "A2": "a2", "A3": "a3", "A4": "a4", "A5": "a5", "A6": "a6",
    # Framework Adaptation: the FOUR retained survey items (Table A-0) -> fa1..fa4; fa5 = expert composite
    "B1": "fa1", "B2": "fa2", "B3": "fa3", "B4": "fa4",          # <-- EDIT to the four retained B items
    # Resource Reorganization survey items (3 per dimension)
    "C1": "drr1", "C2": "drr2", "C3": "drr3",
    "C4": "trr1", "C5": "trr2", "C6": "trr3",
    "C7": "orr1", "C8": "orr2", "C9": "orr3",
    "C10": "crr1", "C11": "crr2", "C12": "crr3",
    # Data Commercialization survey items (3 retained per dimension; EDIT the dropped EV / DIV item)
    "D1": "ev1", "D2": "ev2", "D3": "ev3",
    "D5": "psv1", "D6": "psv2", "D7": "psv3",
    "D8": "cv1", "D9": "cv2", "D10": "cv3",
    "D11": "div1", "D12": "div2", "D13": "div3",
    # marker variable (job satisfaction, 3 MSQ items)
    "JS1": "js1", "JS2": "js2", "JS3": "js3",
}

# --- 2. archival unit-year file: raw column name -> codebook name -----------------------------
ARCHIVAL_COLUMNS = {
    "business_unit": "unit_id", "year": "year",
    "et": "et", "da": "da",
    "pilots": "pr1_pilots", "subsidy_cny": "subsidy_cny", "tax_incentive": "pr1_tax_incentive",
    "n_requirements": "pr2_n_requirements", "penalty_severity": "pr2_penalty_severity", "free_data_scope": "pr2_free_data_scope",
    "employees": "employees", "founding_year": "founding_year", "n_competitors": "ctrl_n_competitors", "gdp_pc": "gdp_pc",
    "adj_province_support": "iv1_adj_province_support", "obs_station_density": "iv2_obs_station_density",
    # archival composite indicators entering the mixed measurement model
    "fa_expert_composite": "fa5",
    "drr_archival": "drr4", "trr_archival": "trr4", "orr_archival": "orr4", "crr_archival": "crr4",
}


def z(s): return (s - s.mean()) / s.std(ddof=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--survey", required=True); ap.add_argument("--sheet", default=0)
    ap.add_argument("--archival", required=True); ap.add_argument("--archival-sheet", default=0)
    a = ap.parse_args()
    sv = pd.read_excel(a.survey, sheet_name=a.sheet).rename(columns=SURVEY_COLUMNS)
    ar = pd.read_excel(a.archival, sheet_name=a.archival_sheet).rename(columns=ARCHIVAL_COLUMNS)
    sv["wave"] = sv["year"] - sv["year"].min() + 1
    # derived archival variables (Section 4.5.1, 4.5.6)
    if "subsidy_cny" in ar: ar["pr1_log_subsidy"] = np.log1p(ar["subsidy_cny"])
    if "employees" in ar: ar["ctrl_log_employees"] = np.log1p(ar["employees"])
    if "founding_year" in ar: ar["ctrl_firm_age"] = ar["year"] - ar["founding_year"]
    if "gdp_pc" in ar: ar["ctrl_log_gdp_pc"] = np.log1p(ar["gdp_pc"])
    ar["PR1"] = z(z(ar["pr1_pilots"]) + z(ar["pr1_log_subsidy"]) + z(ar["pr1_tax_incentive"]))
    ar["PR2"] = z(z(ar["pr2_n_requirements"]) + z(ar["pr2_penalty_severity"]) + z(ar["pr2_free_data_scope"]))
    ar["PR"] = z(ar["PR1"] - ar["PR2"]); ar["PR1xPR2"] = ar["PR1"] * ar["PR2"]
    df = sv.merge(ar, on=["unit_id", "year"], how="left", validate="many_to_one")
    # composites used by the Python cross-checks
    blocks = {"FA": [f"fa{i}" for i in range(1, 6)], "DRR": [f"drr{i}" for i in range(1, 5)], "TRR": [f"trr{i}" for i in range(1, 5)],
              "ORR": [f"orr{i}" for i in range(1, 5)], "CRR": [f"crr{i}" for i in range(1, 5)], "EV": ["ev1", "ev2", "ev3"],
              "PSV": ["psv1", "psv2", "psv3"], "CV": ["cv1", "cv2", "cv3"], "DIV": ["div1", "div2", "div3"]}
    for k, cols in blocks.items():
        have = [c for c in cols if c in df]; df[k] = df[have].mean(axis=1)
    df["RR"] = df[["DRR", "TRR", "ORR", "CRR"]].mean(axis=1); df["DC"] = df[["EV", "PSV", "CV", "DIV"]].mean(axis=1)
    df["PR1_perc"] = df[["a1", "a2", "a3"]].mean(axis=1); df["PR2_perc"] = df[["a4", "a5", "a6"]].mean(axis=1)
    if {"js1", "js2", "js3"} <= set(df): df["JS_marker"] = df[["js1", "js2", "js3"]].mean(axis=1)
    OUT.parent.mkdir(parents=True, exist_ok=True); df.to_csv(OUT, index=False)
    missing = [c for c in ["employee_id", "wave", "unit_id", "et", "da", "PR1", "PR2"] + sum(blocks.values(), []) if c not in df]
    print(f"Wrote {OUT}  ({len(df)} rows, {df.employee_id.nunique()} employees, {df.unit_id.nunique()} units)")
    print("Missing codebook columns:", missing or "none")


if __name__ == "__main__":
    main()
