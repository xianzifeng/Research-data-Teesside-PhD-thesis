# Codebook — `synthetic_panel_person_wave.csv` (and the real panel, same layout)

One row = one employee in one annual survey wave (person-wave). 1,800 rows = 300 employees × 6 waves (2020–2025). All Likert items are integers 1–5 (1 = strongly disagree, 5 = strongly agree); **no item is reverse-coded** (thesis Appendix A). Composite scores are item means. Variables marked *archival* are measured at the business-unit-year level and merged onto every person-wave of that unit and year (54 unit-years; also provided separately in `synthetic_unit_year.csv`).

## Identifiers and design variables

| Variable | Type | Description | Thesis |
|---|---|---|---|
| `employee_id` | string | Pseudonymous respondent id (E001–E300); cluster variable for robust SEs and bootstrap | §4.4.4, §4.9.3 |
| `unit_id`, `unit_label` | int / string | Business unit 1–9 (SUB-A … SUB-I; generic labels, Table 4-2) | Table 4-2 |
| `wave`, `year` | int | Survey wave 1–6 = 2020–2025 | §4.4.1 |
| `role` | string | managerial (180) / technical (120) | §5.1 |
| `tenure_2020` | int | Years with the enterprise at Wave 1 (synthetic only) | — |
| `et` | 0/1 | Enterprise type proxy: 1 = market-facing subsidiary, 0 = internal business unit | §4.5.5, §5.5.2 |
| `da` | 0/1 | Data attribute: 1 = value-added / processed data unit, 0 = basic public data dissemination | §4.5.5, §5.5.2 |
| `region` | string | East / Central / West (Table D-3) | Appendix D.3 |

## Measurement items (Table 4-4 codes ↔ Appendix A codes)

| Construct (dimension) | Items in this file | Table 4-4 codes | Appendix A codes |
|---|---|---|---|
| Policy support perception (PR1) | `a1 a2 a3` | — (combined with archival PR1, §4.4.1 Stage 4) | A1–A3 |
| Compliance constraint perception (PR2; high = heavier constraints) | `a4 a5 a6` | — (combined with archival PR2) | A4–A6 |
| Framework Adaptation (FA) | `fa1 … fa5` | FA1–FA5 | B-block |
| Data resource reorganization (DRR) | `drr1 … drr4` | DRR1–DRR4 | C-block |
| Technology resource reorganization (TRR) | `trr1 … trr4` | TRR1–TRR4 | C-block |
| Organizational resource reorganization (ORR) | `orr1 … orr4` | ORR1–ORR4 | C-block |
| Compliance resource reorganization (CRR) | `crr1 … crr4` | CRR1–CRR4 | C-block |
| Economic value (EV) | `ev1 ev2 ev3` | EV1–EV3 | D-block |
| Public service value (PSV) | `psv1 psv2 psv3` | PSV1–PSV3 | D-block |
| Compliance value (CV) | `cv1 cv2 cv3` | CV1–CV3 | D-block |
| Data integration value (DIV) | `div1 div2 div3` | DIV1–DIV3 | D-block |
| Job satisfaction marker (3 MSQ-type items) | `js1 js2 js3` | — | §5.4 marker variable |

> The exact mapping of the 33 retained Table 4-4 codes onto the 40 administered Appendix A items (and the seven items dropped after the Wave-1 EFA) is to be documented by the author in Appendix A (revision note C379 / C747). Until then the block-level correspondence above applies.

## Composite scores (item means)

`FA`, `DRR`, `TRR`, `ORR`, `CRR`, `RR` (= mean of the four RR dimensions), `EV`, `PSV`, `CV`, `DIV`, `DC` (= mean of the four DC dimensions), `PR1_perc`, `PR2_perc`, `JS_marker`.
Standardised versions `zFA`, `zRR`, `zDC`, `zPR1`, `zPR2`, `zPR1xPR2` are added by `01_data_preprocessing.py` to `outputs/analysis_panel.csv`.

## Archival Policy Regulation composites (Section 4.5.1) — unit-year level

| Variable | Description |
|---|---|
| `pr1_pilots` | Number of government pilot projects the unit participated in (count) |
| `pr1_log_subsidy` | log(1 + policy subsidy amount) |
| `pr1_tax_incentive` | Tax incentive intensity, 1 (none) – 5 (full exemption) |
| `pr2_n_requirements` | Number of mandatory compliance requirements (count) |
| `pr2_penalty_severity` | Severity of violation penalties, 1 (warning) – 5 (revocation) |
| `pr2_free_data_scope` | Scope of data that must be provided free, 1 (basic only) – 5 (all observation data) |
| `PR1` | Standardised policy support index (z of the sum of z-scored PR1 indicators) |
| `PR2` | Standardised compliance constraint index |
| `PR` | Standardised net policy regulation score = z(PR1 − PR2) |
| `PR1xPR2` | Product term for H10 |

## Controls (Section 4.5.6) — unit-year level

`ctrl_log_employees` (log 1+ employees), `ctrl_firm_age` (years since establishment), `ctrl_n_competitors` (local competitors), `ctrl_log_gdp_pc` (log regional GDP per capita).

## Instruments (Section 4.6.1 item 4) — unit-year level

`iv1_adj_province_support` = adjacent-province policy support intensity; `iv2_obs_station_density` = regional observation-station density.

## Archival outcome indicators (Section 4.5.4; mixed measurement model) — unit-year level

`arch_dc_revenue_growth_pct` = data-commercialization revenue growth (%); `arch_dc_direct_sales_share_pct` = direct data sales / data-commercialization revenue (%).

## Mplus names (`mplus/panel.dat`)

Same variables, truncated to 8 characters: `emp` (numeric employee id), `unit`, `wave`, `et`, `da`, items as above, `js1–js3`, `pr1`, `pr2`, `pr1xpr2`, `iv1`, `iv2`, `lnemp`, `age`, `ncomp`, `lngdp`. Missing code 9999. The authoritative order is written to `mplus/mplus_variable_order.txt` on every export.
