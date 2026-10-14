"""
Shared definitions for the version 4 analyses.

Each agency is judged by its own intensity categories, at three thresholds:

    all   the definition used in multiverse.py: a JMA grade of 2-5 or 9, or any
          positive wind from JTWC, CMA or HKO
    TD    tropical depression or stronger: JMA grade 2-5 or 9, USA_SSHS >= -1,
          CMA_CAT 1-6, HKO_CAT TD or stronger
    TS    tropical storm or stronger: JMA grade 3-5 or 9, USA_SSHS >= 0,
          CMA_CAT 2-6, HKO_CAT TS or stronger

A storm meets a threshold for an agency if that agency's category meets it at
any fix inside the domain, as in multiverse.py. CMA_TSwind is a sensitivity
flag that classifies CMA by wind (CMA_WIND >= 34 kt, extratropical stage
excluded) instead of by category.

The rebuilt grading requirements use the three basin-wide agencies: JMA, JTWC
and CMA alone, JMA or JTWC, any of the three, and all three. HKO enters the
archive in 1961 and is treated separately.

fit_many applies the test in multiverse.fit to every row of a matrix at once;
check_v4.py confirms it reproduces multiverse.fit verdict for verdict.
"""
import numpy as np
import pandas as pd
from scipy import stats

from multiverse import FIRST, LAST

YEARS = list(range(FIRST, LAST + 1))
AGENCIES = ("JMA", "JTWC", "CMA", "HKO")
THRESHOLDS = ("all", "TD", "TS")
REQUIREMENTS = ("JMA", "JTWC", "CMA", "JMA|JTWC", "any of 3", "all 3")
FLAGS = [f"{a}_{t}" for a in AGENCIES for t in THRESHOLDS] + ["CMA_TSwind"]


def _num(d, c):
    return pd.to_numeric(d[c], errors="coerce")


def masks(d):
    """Row-level flags for every agency and threshold."""
    g = _num(d, "TOKYO_GRADE" if "TOKYO_GRADE" in d.columns else "TOK_GRADE")
    usw, ss = _num(d, "USA_WIND").fillna(0), _num(d, "USA_SSHS")
    cw, cc = _num(d, "CMA_WIND").fillna(0), _num(d, "CMA_CAT")
    hw, hc = _num(d, "HKO_WIND").fillna(0), d["HKO_CAT"].astype(str).str.strip()
    m = {
        "JMA_all": g.isin([2, 3, 4, 5, 9]),
        "JMA_TD": g.isin([2, 3, 4, 5, 9]),
        "JMA_TS": g.isin([3, 4, 5, 9]),
        "JTWC_all": usw > 0,
        "JTWC_TD": ss >= -1,
        "JTWC_TS": ss >= 0,
        "CMA_all": cw > 0,
        "CMA_TD": cc.isin([1, 2, 3, 4, 5, 6]),
        "CMA_TS": cc.isin([2, 3, 4, 5, 6]),
        "HKO_all": hw > 0,
        "HKO_TD": hc.isin(["TD", "TS", "STS", "T", "ST", "SuperT"]),
        "HKO_TS": hc.isin(["TS", "STS", "T", "ST", "SuperT"]),
        "CMA_TSwind": (cw >= 34) & (cc != 9),
    }
    return {k: v.fillna(False).astype(bool) for k, v in m.items()}


def storm_flags(sub):
    """One row per SID: its season and whether each flag holds at any fix."""
    mk = masks(sub)
    f = pd.DataFrame({k: mk[k].values for k in FLAGS}, index=sub.index)
    f["SID"] = sub.SID.values
    f["SEASON"] = sub.SEASON.values
    g = f.groupby("SID")
    se = g.SEASON.agg(["min", "max"])
    if not (se["min"] == se["max"]).all():
        raise ValueError("SEASON varies within a storm")
    out = g[FLAGS].any()
    out.insert(0, "SEASON", se["min"])
    return out


def requirement(fl, req, t):
    J, T, C = fl[f"JMA_{t}"], fl[f"JTWC_{t}"], fl[f"CMA_{t}"]
    return {"JMA": J, "JTWC": T, "CMA": C, "JMA|JTWC": J | T,
            "any of 3": J | T | C, "all 3": J & T & C}[req]


def any4(fl, t):
    return fl[f"JMA_{t}"] | fl[f"JTWC_{t}"] | fl[f"CMA_{t}"] | fl[f"HKO_{t}"]


def annual(fl, sel):
    """Storms per season, 1884 to 2023, for the storms selected by sel."""
    return fl[sel].groupby("SEASON").size().reindex(YEARS).fillna(0).astype(int)


def load_flags(path, domain, rule="main track"):
    f = pd.read_csv(path)
    f = f[(f.domain == domain) & (f.rule == rule)].set_index("SID")
    for c in FLAGS:
        f[c] = f[c].astype(bool)
    return f


def fit_many(V, x):
    """The test in multiverse.fit, applied to every row of V (series x years).
    Returns slope, half-width of the no-trend band, and verdict (-1, 0, +1)."""
    V = np.asarray(V, float)
    x = np.asarray(x, float)
    n = V.shape[1]
    xc = x - x.mean()
    den = (xc ** 2).sum()
    slope = V @ xc / den
    R = V - V.mean(1, keepdims=True) - slope[:, None] * xc[None, :]
    a, c = R[:, :-1], R[:, 1:]
    am, cm = a.mean(1, keepdims=True), c.mean(1, keepdims=True)
    num = ((a - am) * (c - cm)).sum(1)
    dd = np.sqrt(((a - am) ** 2).sum(1) * ((c - cm) ** 2).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r1 = np.clip(num / dd, 0.0, 0.99)          # NaN stays NaN -> verdict 0
        ne = np.maximum(n * (1 - r1) / (1 + r1), 4.0)
        se = np.sqrt((R ** 2).sum(1) / (n - 2) / den) * np.sqrt(n / ne)
        band = stats.t.ppf(0.975, np.maximum(ne - 2, 1)) * se
        verdict = np.where(slope < -band, -1, np.where(slope > band, 1, 0))
    return slope, band, verdict
