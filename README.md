# Replication materials — *Institutional Constraints and Public Data Commercialization: Evidence from a Chinese Meteorological Media Enterprise*

Doctoral thesis, Teesside University International Business School (module OBE5008-N).
Author: Bai Yu · Supervisor: Murphy · Submitted October 2026.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23448040.svg)](https://doi.org/10.5281/zenodo.23448040)  Archived release 1.0.1 (18 October 2026).

This repository contains everything needed to re-run the quantitative analyses of Chapters 4–5 and the MMDV valuation worksheet of Appendix E:

* the **Python code** behind Appendix C (preprocessing, measurement checks, mediation bootstrap, IV estimation, moderation, robustness, MMDV);
* the **Mplus 8.8 input files** for the CFA, structural SEM, bootstrap mediation, PR1 × PR2 interaction and multi-group (ET / DA) models;
* a **synthetic dataset** with exactly the structure of the Huafeng panel (300 employees × 6 waves, 9 business units, 33 + 6 Likert items, archival PR1/PR2 composites, instruments, controls, marker items) so that every script runs end-to-end;
* `requirements.txt`, a codebook and this documentation.

> **Data availability.** The real person-wave panel and the business-unit financial indicators are **not** deposited: they are covered by the organisational access approval and data-sharing agreements described in thesis Appendix H.2. They are classified as confidential internal enterprise data and are **not provided outside the research team in any form**, including de-identified extracts (thesis §4.8.4, Appendix H.4). The synthetic file reproduces the variable structure and the target moments reported in Tables 5-1 to 5-3; it is **not** real data and point estimates obtained from it will differ from the thesis.

## Repository structure

```
Research-data-Teesside-PhD-thesis/
├── README.md                      ← this file (EN); 中文说明见 docs/STRUCTURE_zh.md
├── LICENSE                        ← MIT (code) / CC-BY-4.0 (synthetic data & docs)
├── CITATION.cff                   ← citation metadata, read automatically by GitHub and Zenodo
├── requirements.txt               ← pinned Python environment (Python 3.11)
├── .gitignore
│
├── code/                          ← Python scripts (thesis Appendix C), run in numerical order
│   ├── 00_generate_synthetic_data.py   builds data/synthetic/* from the thesis parameters
│   ├── 01_data_preprocessing.py        Appendix C.1 — quality checks, VIF, outliers (§4.4.3, §5.4)
│   ├── 02_measurement_checks.py        α / CR / AVE / HTMT / EFA / KMO / ICC / Harman (§4.5.7, §4.7.1, §5.4, Table D-1)
│   ├── 03_structural_paths_crosscheck.py  Appendix C.2 cross-check of Table 5-3 (cluster-robust OLS, two-way FE)
│   ├── 04_bootstrap_mediation.py       Appendix C.3 — statsmodels Mediation + employee-cluster bootstrap (Tables 5-4, 5-16)
│   ├── 05_iv_estimation.py             Appendix C.4 — 2SLS with IV1/IV2 (Table 5-3 Panel B)
│   ├── 06_moderation_interaction.py    H8–H10: PR1×PR2, simple slopes, Johnson–Neyman, ET/DA groups (Tables 5-6, 5-7, 5-13)
│   ├── 07_robustness_checks.py         Appendix C.5 — leave-one-out, placebo, winsorising, regions, wild-cluster bootstrap
│   ├── 08_mmdv_worksheet.py            Appendix E — CBV / MBV / IBV, composite, sensitivity (Tables 5-14, 5-15)
│   ├── 09_prepare_mplus_data.py        exports mplus/panel.dat for the .inp files
│   ├── 10_sem_latent_semopy.py         open-source latent SEM equivalent of mplus/01-05 (semopy + employee-cluster bootstrap)
│   ├── 11_prepare_real_panel.py        converts the author's raw files into the codebook layout (run locally only)
│   ├── run_all.py                      runs 00–09 and writes outputs/run_log.txt
│   └── utils.py                        shared helpers (α, CR/AVE, HTMT, ICC, KMO, cluster bootstrap, BC intervals)
│
├── mplus/                         ← Mplus 8.8 syntax (see mplus/README.md)
│   ├── 01_cfa_measurement.inp
│   ├── 02_sem_structural_mlr.inp
│   ├── 03_mediation_bootstrap.inp
│   ├── 04_interaction_PR1xPR2.inp
│   ├── 05_multigroup_ET_DA.inp
│   ├── panel.dat                       generated; header-less data for Mplus
│   └── mplus_variable_order.txt        generated; NAMES statement that matches panel.dat
│
├── data/
│   ├── README.md                       data statement
│   ├── synthetic/
│   │   ├── synthetic_panel_person_wave.csv   1,800 person-waves × 90 variables (generated)
│   │   ├── synthetic_unit_year.csv           54 business-unit-years of archival variables (generated)
│   │   ├── generation_parameters.json        every target parameter with its thesis source
│   │   └── generation_log.txt                achieved vs. target moments
│   └── real/
│       └── README.md                   how to request the confidential panel (nothing else is stored here)
│
├── docs/
│   ├── CODEBOOK.md                     variable-by-variable codebook (items ↔ Table 4-4 / Appendix A)
│   └── STRUCTURE_zh.md                 仓库结构与使用说明（中文）
│
├── notebooks/
│   └── 01_replication_walkthrough.ipynb   runs the pipeline step by step with commentary
│
└── outputs/                       ← created by the scripts (csv tables + run_log.txt); safe to delete
```

## Quick start

```bash
git clone https://github.com/xianzifeng/Research-data-Teesside-PhD-thesis.git
cd Research-data-Teesside-PhD-thesis
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python code/run_all.py                                  # ~3-5 min; bootstrap steps dominate
```

`run_all.py` regenerates the synthetic data, runs every Python check and writes result tables to `outputs/` and a full console log to `outputs/run_log.txt`. To run the Mplus models, open the `.inp` files in `mplus/` after step 09 has written `panel.dat`.

If the confidential panel is available (`data/real/huafeng_panel_person_wave.csv`, same column names as the synthetic file), the scripts pick it up automatically and skip the synthetic generation.

## Mapping between the thesis and this repository

| Thesis element | Where it is reproduced |
|---|---|
| Table 4-4 / 5-2 (loadings, α, CR, AVE) | Mplus `01_cfa_measurement.inp`; Python approximation in `02_measurement_checks.py` |
| Table D-1 (Wave-1 EFA), KMO/Bartlett, parallel analysis | `02_measurement_checks.py` |
| Table 5-1 (descriptives, correlations) | `generation_log.txt`; `02_measurement_checks.py` |
| Table 5-3 (structural paths) | Mplus `02_sem_structural_mlr.inp`; composite cross-check `03_structural_paths_crosscheck.py` |
| Table 5-3 Panel B (IV) | `05_iv_estimation.py` |
| Tables 5-4, 5-16 (mediation, bootstrap CIs) | Mplus `03_mediation_bootstrap.inp`; `04_bootstrap_mediation.py` |
| Tables 5-6, 5-7 (ET / DA groups), Table 5-13 (simple slopes), Johnson–Neyman | Mplus `04_…`, `05_…`; `06_moderation_interaction.py` |
| Tables 5-5, 5-8 to 5-12, D-3, D-4 (robustness) | `07_robustness_checks.py` |
| Section 5.4 (VIF, Cook's D, Mahalanobis, Harman, marker variable) | `01_data_preprocessing.py`, `02_measurement_checks.py` |
| Section 4.9.3 (ICC, effective N) | `02_measurement_checks.py` |
| Appendix E / Tables 5-14, 5-15 (MMDV) | `08_mmdv_worksheet.py` (reproduces the worksheet figures exactly) |

## Running the latent-variable models without Mplus

`code/10_sem_latent_semopy.py` fits the same five models as the Mplus inputs with **semopy** (`pip install semopy`) and replaces TYPE = COMPLEX / BOOTSTRAP with an employee-cluster bootstrap (whole employees resampled), reporting standardised paths, fit indices (χ², CFI, TLI, RMSEA, SRMR), CR/AVE, indirect effects with percentile and bias-corrected intervals, the PR1 × PR2 simple slopes and the ET/DA group comparison. On the real panel: `python code/11_prepare_real_panel.py --survey … --archival …` then `python code/10_sem_latent_semopy.py --data data/real/huafeng_panel_person_wave.csv --boot 5000 --out outputs/real_results`. Only summary tables are written; no row-level data leave the computer.

## Software versions

The thesis analyses were run with Python 3.11 — pandas 2.0.3, numpy 1.24.3, statsmodels 0.14.0, linearmodels 4.31, scipy 1.11.1, scikit-learn 1.3.0 — and Mplus 8.8. `requirements.txt` pins these versions; the scripts also run on current releases (tested with pandas 3.0 / statsmodels 0.15 / linearmodels 7.1).

## Citing

See `CITATION.cff`. Suggested form:

> Bai, Y. (2026) *Replication materials for 'Institutional constraints and public data commercialization: evidence from a Chinese meteorological media enterprise'* (Version 1.0.1) [Software and synthetic data set]. Zenodo. https://doi.org/10.5281/zenodo.23448040

## Licence

Code: MIT. Synthetic data, codebook and documentation: CC BY 4.0. No real respondent or enterprise data are included.
