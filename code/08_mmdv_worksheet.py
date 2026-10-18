"""
08_mmdv_worksheet.py  (thesis Appendix E; Section 5.5.5; Tables 5-14, 5-15)
Reproduces the Multi-Method Meteorological Data Valuation (MMDV) worksheet for
the WMO/NOAA proxy portfolio step by step, including the sensitivity scenarios.
All inputs are the published benchmark figures and expert-elicited parameters
listed in Appendix E; change them in PARAMS to run alternative scenarios.
This is an illustrative valuation of a public proxy portfolio, not of Huafeng's assets.
"""
import json
import pandas as pd

from utils import OUT

PARAMS = {
    # E.1 cost-based
    "capex_gbon_bn": 2.0, "opex_annual_bn": 1.5, "horizon_years": 3,
    "timeliness_decay": 0.15, "spatial_premium": 1.25, "scenario_multiplier": 1.40,
    # E.2 market-based (annual benefit flows, US$ bn)
    "comparables_bn": {"NOAA CONUS radar network": 1.176, "WMO surface observation improvement (global)": 5.19,
                       "NOAA radar outage mitigation": 0.0291, "US Landsat archive": 1.7},
    "china_gdp_share_mbv": 0.15, "quality_adjustment": 0.90, "domestic_premium": 1.10,
    # E.3 income-based
    "wmo_total_benefit_bn": 162.0,
    "sector_benefit_bn": {"disaster management": 45, "agriculture": 38, "energy": 32, "transportation": 28, "water supply": 12, "construction": 7},
    "sector_capture_rate": {"disaster management": 0.03, "agriculture": 0.08, "energy": 0.10, "transportation": 0.06, "water supply": 0.05, "construction": 0.04},
    "china_gdp_share_ibv": 0.18, "risk_discount": 0.25,
    # E.4 weights
    "weights": {"CBV": 0.30, "MBV": 0.35, "IBV": 0.35},
}


def mmdv(p, exclude_global=False):
    out = {}
    cost = p["capex_gbon_bn"] + p["opex_annual_bn"] * p["horizon_years"]
    out["replacement_cost"] = cost
    out["CBV"] = cost * (1 - p["timeliness_decay"]) * p["spatial_premium"] * p["scenario_multiplier"]
    comps = {k: v for k, v in p["comparables_bn"].items() if not (exclude_global and "global" in k)}
    avg = sum(comps.values()) / len(comps)
    out["MBV_avg_annual"] = avg
    out["MBV"] = avg * p["china_gdp_share_mbv"] * p["quality_adjustment"] * p["domestic_premium"] * p["horizon_years"]
    num = sum(p["sector_benefit_bn"][s] * p["sector_capture_rate"][s] for s in p["sector_benefit_bn"])
    cap = num / sum(p["sector_benefit_bn"].values())
    out["capture_rate"] = cap
    annual = p["wmo_total_benefit_bn"] * p["china_gdp_share_ibv"] * cap
    out["IBV_annual"] = annual
    out["IBV"] = annual * p["horizon_years"] * (1 - p["risk_discount"])
    out["composite_equal"] = (out["CBV"] + out["MBV"] + out["IBV"]) / 3
    w = p["weights"]
    out["composite_weighted"] = out["CBV"] * w["CBV"] + out["MBV"] * w["MBV"] + out["IBV"] * w["IBV"]
    out["benefit_cost_ratio_MBV_over_cost"] = out["MBV"] / cost
    return out


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    base = mmdv(PARAMS)
    print("MMDV base case (US$ billion, 3-year horizon)")
    for k, v in base.items():
        print(f"  {k:32s} {v:8.3f}")
    print("\nThesis reference: CBV 9.67 | MBV 0.90 | IBV 4.11 | equal 4.89 | weighted 4.65 | BCR 0.14")
    scen = {
        "Conservative": dict(timeliness_decay=0.25, spatial_premium=1.10, scenario_multiplier=1.20, risk_discount=0.35),
        "Base": {},
        "Optimistic": dict(timeliness_decay=0.10, spatial_premium=1.40, scenario_multiplier=1.60, risk_discount=0.15),
    }
    rows = []
    for name, over in scen.items():
        p = json.loads(json.dumps(PARAMS)); p.update(over)
        r = mmdv(p)
        rows.append({"scenario": name, **{k: p[k] for k in ["timeliness_decay", "spatial_premium", "scenario_multiplier", "risk_discount"]},
                     "CBV": r["CBV"], "MBV": r["MBV"], "IBV": r["IBV"], "weighted_composite": r["composite_weighted"]})
    tab = pd.DataFrame(rows)
    print("\nSensitivity (Table 5-14 / E.4 Step 3):")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    print(f"\nMBV excluding the global aggregate: {mmdv(PARAMS, exclude_global=True)['MBV']:.3f} bn  (thesis 0.431)")
    tab.to_csv(OUT / "08_mmdv_sensitivity.csv", index=False)
