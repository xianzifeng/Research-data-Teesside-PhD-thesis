"""
09_prepare_mplus_data.py
Exports the person-wave panel to the header-less, space-delimited .dat file
read by every input in /mplus (Mplus 8.8 does not read column headers).
The variable order written here MUST match the NAMES statement in the .inp files;
the order is also written to mplus/mplus_variable_order.txt for checking.
Missing values are coded 9999 (MISSING ARE ALL (9999) in the .inp files).
Employee id is exported as an integer (Mplus CLUSTER requires numeric ids).
"""
import pandas as pd
from utils import load_panel, BLOCKS, CONTROLS, ROOT

df = load_panel()
df["emp_num"] = pd.factorize(df.employee_id)[0] + 1
cols = (["emp_num", "unit_id", "wave", "et", "da"]
        + sum([BLOCKS[b] for b in ["PR1_perc", "PR2_perc", "FA", "DRR", "TRR", "ORR", "CRR", "EV", "PSV", "CV", "DIV"]], [])
        + ["js1", "js2", "js3", "PR1", "PR2", "PR1xPR2", "iv1_adj_province_support", "iv2_obs_station_density"] + CONTROLS)
mp = df[cols].copy()
mp.columns = [c[:8] for c in mp.columns]  # Mplus variable names: max 8 characters
rename = {"iv1_adj_": "iv1", "iv2_obs_": "iv2", "ctrl_log": "lnemp", "ctrl_fir": "age", "ctrl_n_c": "ncomp"}
# the two ctrl_log_* names collide after truncation -> assign explicitly
mp.columns = ["emp", "unit", "wave", "et", "da"] + cols[5:5 + 39] + ["js1", "js2", "js3", "pr1", "pr2", "pr1xpr2", "iv1", "iv2", "lnemp", "age", "ncomp", "lngdp"]
assert len(set(mp.columns)) == mp.shape[1] and all(len(c) <= 8 for c in mp.columns)
(ROOT / "mplus").mkdir(exist_ok=True)
mp.round(4).fillna(9999).to_csv(ROOT / "mplus" / "panel.dat", sep=" ", header=False, index=False)
(ROOT / "mplus" / "mplus_variable_order.txt").write_text(
    "NAMES = " + " ".join(mp.columns) + ";\n(source: %s; N = %d rows)\n" % (df.attrs["source"], len(mp)))
print(f"Wrote mplus/panel.dat ({len(mp)} rows x {mp.shape[1]} vars) from {df.attrs['source']} data")
print("NAMES =", " ".join(mp.columns))
