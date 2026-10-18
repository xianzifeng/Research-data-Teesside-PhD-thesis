"""Shared helpers for the replication scripts (loading, reliability, validity, clustering)."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_SYNTH = ROOT / "data" / "synthetic" / "synthetic_panel_person_wave.csv"
DATA_REAL = ROOT / "data" / "real" / "huafeng_panel_person_wave.csv"   # not distributed (Appendix H.2)
OUT = ROOT / "outputs"

# Measurement blocks of Table 4-4 (33 retained items) + PR perception items (Appendix A, A1-A6)
BLOCKS = {
    "FA": ["fa1", "fa2", "fa3", "fa4", "fa5"],
    "DRR": ["drr1", "drr2", "drr3", "drr4"],
    "TRR": ["trr1", "trr2", "trr3", "trr4"],
    "ORR": ["orr1", "orr2", "orr3", "orr4"],
    "CRR": ["crr1", "crr2", "crr3", "crr4"],
    "EV": ["ev1", "ev2", "ev3"],
    "PSV": ["psv1", "psv2", "psv3"],
    "CV": ["cv1", "cv2", "cv3"],
    "DIV": ["div1", "div2", "div3"],
    "PR1_perc": ["a1", "a2", "a3"],
    "PR2_perc": ["a4", "a5", "a6"],
}
SECOND_ORDER = {"RR": ["DRR", "TRR", "ORR", "CRR"], "DC": ["EV", "PSV", "CV", "DIV"]}
CONTROLS = ["ctrl_log_employees", "ctrl_firm_age", "ctrl_n_competitors", "ctrl_log_gdp_pc"]


def load_panel(path=None):
    """Load the real panel if present (never distributed), otherwise the synthetic one."""
    if path is None:
        path = DATA_REAL if DATA_REAL.exists() else DATA_SYNTH
    df = pd.read_csv(path)
    df.attrs["source"] = "REAL" if Path(path) == DATA_REAL else "SYNTHETIC"
    return df


def cronbach_alpha(items: pd.DataFrame) -> float:
    k = items.shape[1]
    return k / (k - 1) * (1 - items.var(ddof=1).sum() / items.sum(axis=1).var(ddof=1))


def standardized_loadings(items: pd.DataFrame) -> np.ndarray:
    """One-factor loadings = correlations with the first principal component (close to CFA loadings)."""
    z = (items - items.mean()) / items.std(ddof=1)
    u, s, vt = np.linalg.svd(z.values, full_matrices=False)
    pc1 = u[:, 0] * s[0]
    lam = np.array([np.corrcoef(pc1, z[c])[0, 1] for c in z.columns])
    return np.abs(lam)


def cr_ave(lam: np.ndarray):
    """Composite reliability and average variance extracted from standardised loadings."""
    err = 1 - lam ** 2
    cr = lam.sum() ** 2 / (lam.sum() ** 2 + err.sum())
    return cr, (lam ** 2).mean()


def htmt(df: pd.DataFrame, block_a, block_b) -> float:
    """Heterotrait-monotrait ratio (Henseler, Ringle & Sarstedt, 2015)."""
    R = df[block_a + block_b].corr().abs()
    hetero = R.loc[block_a, block_b].values.mean()
    def mono(b):
        m = R.loc[b, b].values
        return m[np.triu_indices(len(b), 1)].mean()
    return hetero / np.sqrt(mono(block_a) * mono(block_b))


def icc1(y: pd.Series, groups: pd.Series) -> float:
    d = pd.DataFrame({"y": y.values, "g": groups.values})
    m = d.groupby("g").y.agg(["mean", "count"])
    gm = d.y.mean()
    msb = (m["count"] * (m["mean"] - gm) ** 2).sum() / (len(m) - 1)
    msw = ((d.y - d.g.map(m["mean"])) ** 2).sum() / (len(d) - len(m))
    k = m["count"].mean()
    return (msb - msw) / (msb + (k - 1) * msw)


def kmo_bartlett(items: pd.DataFrame):
    """Kaiser-Meyer-Olkin measure and Bartlett's test of sphericity."""
    from scipy import stats
    R = items.corr().values
    inv = np.linalg.pinv(R)
    d = np.sqrt(np.diag(inv))
    P = -inv / np.outer(d, d)           # partial correlations
    np.fill_diagonal(P, 0); R0 = R.copy(); np.fill_diagonal(R0, 0)
    kmo = (R0 ** 2).sum() / ((R0 ** 2).sum() + (P ** 2).sum())
    n, p = items.shape
    chi2 = -(n - 1 - (2 * p + 5) / 6) * np.log(np.linalg.det(R))
    dof = p * (p - 1) / 2
    return kmo, chi2, dof, stats.chi2.sf(chi2, dof)


def cluster_bootstrap_indices(rng, clusters: np.ndarray):
    """Resample whole clusters (employees) with replacement; returns row indices."""
    uniq = np.unique(clusters)
    pick = rng.choice(uniq, size=len(uniq), replace=True)
    pos = {c: np.flatnonzero(clusters == c) for c in uniq}
    return np.concatenate([pos[c] for c in pick])


def bc_interval(boot: np.ndarray, est: float, alpha=0.05):
    """Bias-corrected (BC) bootstrap percentile interval (Efron 1987), as reported by Mplus CINTERVAL(BCBOOTSTRAP)."""
    from scipy import stats
    z0 = stats.norm.ppf((boot < est).mean().clip(1e-6, 1 - 1e-6))
    lo = stats.norm.cdf(2 * z0 + stats.norm.ppf(alpha / 2))
    hi = stats.norm.cdf(2 * z0 + stats.norm.ppf(1 - alpha / 2))
    return np.percentile(boot, 100 * lo), np.percentile(boot, 100 * hi)


def banner(title: str, df: pd.DataFrame):
    print("=" * 78)
    print(f"{title}   [data: {df.attrs.get('source', '?')}, N = {len(df):,}]")
    print("=" * 78)
