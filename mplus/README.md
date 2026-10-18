# Mplus 8.8 input files

| File | Thesis location | Purpose |
|---|---|---|
| `01_cfa_measurement.inp` | §4.4.1 Stage 4, §4.5.7, Tables 4-4 / 5-2 | Second-order CFA of the 33 retained items + A1–A6; TYPE = COMPLEX, MLR; saves factor scores |
| `02_sem_structural_mlr.inp` | Table 5-3, H1–H5 | Structural model PR → FA → RR → DC with controls; MODEL INDIRECT for H6/H7 point estimates |
| `03_mediation_bootstrap.inp` | Table 5-4, Table 5-16, §4.6.1 | 5,000 bootstrap draws, bias-corrected intervals (`CINTERVAL (BCBOOTSTRAP)`); ML without clustering (Mplus restriction, see header notes) |
| `04_interaction_PR1xPR2.inp` | Table 5-13, H10 | PR1 × PR2 interaction on FA (and DC); simple slopes via MODEL CONSTRAINT |
| `05_multigroup_ET_DA.inp` | Tables 5-6 / 5-7, H8 / H9 | Multi-group model by ET (swap to DA); Wald test and indirect-effect difference |
| `panel.dat` | — | Header-less data file written by `code/09_prepare_mplus_data.py` (synthetic unless the real panel is present) |
| `mplus_variable_order.txt` | — | The `NAMES` statement matching `panel.dat` — check it after re-exporting |

Run order: `python code/09_prepare_mplus_data.py` → open each `.inp` in Mplus 8.8 → *Run*.
Output `.out` files are not committed (see `.gitignore`); add them to archive a complete run.

Variable names are truncated to Mplus' 8-character limit:
`lnemp` = log employees, `age` = firm age, `ncomp` = number of competitors, `lngdp` = log regional GDP per capita,
`iv1` = adjacent-province policy support, `iv2` = observation-station density.
