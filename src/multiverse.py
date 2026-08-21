"""
Enumerate every defensible specification of a tropical cyclone count trend
and fit each one.

Three choices enter any such trend, and no convention settles any of them:

    grading req.   which operational agencies must have graded a storm
                   for it to count (6 options, including the raw total
                   that asks for no intensity at all)
    counting rule  which track entries constitute a storm (4 options)
    start year     the first year of the fitted window (50 options,
                   1951 to 2000; every window ends in 2023)

6 x 4 x 50 = 1,200 specifications per domain, fitted for the western North
Pacific basin and again for the PAR polygon, so 2,400 in total. Each slope is
tested against zero using a standard error corrected for serial dependence.

Usage
    python src/multiverse.py ibtracs.WP.list.v04r01.csv par_clipped.csv \
                             data/multiverse.csv
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

FIRST, LAST = 1884, 2023
STARTS = range(1951, 2001)
AGENCIES = ["raw total", "all four", "JMA+JTWC+CMA", "JMA", "JMA+JTWC", "JTWC"]
YEARS = list(range(FIRST, LAST + 1))


def num(df, col):
    return pd.to_numeric(df[col], errors="coerce").fillna(0)


def counting_rules(d):
    """The four ways of deciding which track entries constitute a storm."""
    main = d[d.TRACK_TYPE == "main"]
    fixes = main.groupby("SID").size()
    # IBTrACS NATURE: TS is tropical, DS is disturbance, so a tropical stage is TS alone.
    trop = main.groupby("SID").NATURE.agg(lambda v: "TS" in set(v))
    return {"all SIDs": d,
            "main track": main,
            ">=2 fixes": main[main.SID.isin(fixes[fixes >= 2].index)],
            "tropical nature": main[main.SID.isin(trop[trop].index)]}


def agency_mask(d, key):
    """JMA grades: 2 TD, 3 TS, 4 STS, 5 TY, 9 TC of TS intensity or higher.
    6 is extratropical and 7 is a position marker, so neither is a grade."""
    gcol = "TOKYO_GRADE" if "TOKYO_GRADE" in d.columns else "TOK_GRADE"
    jma = num(d, gcol).isin([2, 3, 4, 5, 9])
    jtwc, cma, hko = num(d, "USA_WIND") > 0, num(d, "CMA_WIND") > 0, num(d, "HKO_WIND") > 0
    return {"raw total": pd.Series(True, index=d.index),
            "JTWC": jtwc, "JMA": jma, "JMA+JTWC": jma | jtwc,
            "JMA+JTWC+CMA": jma | jtwc | cma,
            "all four": jma | jtwc | cma | hko}[key]


def annual(d, mask):
    return (d[mask].groupby("SEASON").SID.nunique()
            .reindex(YEARS).fillna(0).astype(int))


def fit(s, key=None):
    """OLS slope, its standard error corrected for serial dependence, and a
    verdict at the 5% level.

    The correction follows the standard effective-sample-size approach: the
    lag-one autocorrelation of the residuals about the fitted line reduces the
    number of independent observations to n(1-r)/(1+r), and the standard error
    and degrees of freedom are scaled accordingly.

    This replaces a circular block bootstrap used in an earlier version of this
    work. That test was measured, on synthetic series with no trend, to reject
    at 10 to 15 percent against a nominal 5 percent, because residuals that are
    orthogonal to the year index by construction lose variance when resampled in
    long blocks. The estimator below was measured at 4.3 percent for serially
    independent counts, 6.6 percent at a lag-one autocorrelation of 0.33, and
    9.1 percent at 0.50. See src/size.py.
    """
    v = np.asarray(s, float)
    x = np.asarray(s.index, float)
    xc = x - x.mean()
    den = (xc ** 2).sum()
    n = len(v)
    slope = (v @ xc) / den
    resid = v - np.poly1d(np.polyfit(x, v, 1))(x)
    r1 = np.corrcoef(resid[:-1], resid[1:])[0, 1]
    r1 = min(max(r1, 0.0), 0.99)
    n_eff = max(n * (1 - r1) / (1 + r1), 4.0)
    se = np.sqrt((resid @ resid) / (n - 2) / den) * np.sqrt(n / n_eff)
    crit = stats.t.ppf(0.975, max(n_eff - 2, 1))
    lo, hi = -crit * se, crit * se
    return slope, lo, hi, ("+" if slope > hi else "-" if slope < lo else "0")


def main():
    wnp = sys.argv[1] if len(sys.argv) > 1 else "ibtracs.WP.list.v04r01.csv"
    par = sys.argv[2] if len(sys.argv) > 2 else "par_clipped.csv"
    dst = sys.argv[3] if len(sys.argv) > 3 else "data/multiverse.csv"

    basin = pd.read_csv(wnp, low_memory=False, skiprows=[1])
    clip = pd.read_csv(par, low_memory=False)
    rows = []
    for domain, d in [("WNP basin", basin), ("PAR", clip)]:
        d = d[(d.SEASON >= FIRST) & (d.SEASON <= LAST)].copy()
        for rule, sub in counting_rules(d).items():
            for a in AGENCIES:
                series = annual(sub, agency_mask(sub, a))
                for y in STARTS:
                    slope, lo, hi, verdict = fit(series.loc[y:LAST])
                    rows.append(dict(domain=domain, rule=rule, agency=a, start=y,
                                     slope=slope, lo=lo, hi=hi, verdict=verdict))
    out = pd.DataFrame(rows)
    out.to_csv(dst, index=False)
    for domain, g in out.groupby("domain", sort=False):
        dec, inc, n = (g.verdict == "-").sum(), (g.verdict == "+").sum(), len(g)
        print(f"{domain:10s} n={n}  decline {dec} ({100*dec/n:.0f}%)  "
              f"no trend {n-dec-inc} ({100*(n-dec-inc)/n:.0f}%)  increase {inc}")
    print(f"written: {dst}")


if __name__ == "__main__":
    main()
