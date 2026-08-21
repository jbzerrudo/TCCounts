"""
Enumerate every defensible specification of a tropical cyclone count trend
and fit each one.

Three choices enter any such trend, and no convention settles any of them:

    agency set     which operational agencies must have graded a storm
                   for it to count (6 options, including the raw total
                   that asks for no intensity at all)
    counting rule  which track entries constitute a storm (4 options)
    start year     the first year of the fitted window (50 options,
                   1951 to 2000; every window ends in 2023)

6 x 4 x 50 = 1,200 specifications per domain, fitted for the western North
Pacific basin and again for the PAR polygon, so 2,400 in total. Each is
tested against its own circular block bootstrap null.

Usage
    python src/multiverse.py ibtracs.WP.list.v04r01.csv par_clipped.csv \
                             data/multiverse.csv
"""
import sys
import numpy as np
import pandas as pd

FIRST, LAST = 1884, 2023
STARTS = range(1951, 2001)
SEED, NB, BLOCK = 41, 3000, 10
AGENCIES = ["raw total", "all four", "JMA+JTWC+CMA", "JMA", "JMA+JTWC", "JTWC"]
YEARS = list(range(FIRST, LAST + 1))
rng = np.random.default_rng(SEED)


def num(df, col):
    return pd.to_numeric(df[col], errors="coerce").fillna(0)


def counting_rules(d):
    """The four ways of deciding which track entries constitute a storm."""
    main = d[d.TRACK_TYPE == "main"]
    fixes = main.groupby("SID").size()
    trop = main.groupby("SID").NATURE.agg(lambda v: bool(set(v) & {"TS", "DS"}))
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


def fit(s, L=BLOCK, nb=NB):
    """OLS slope, and the 2.5th to 97.5th percentile band of slopes from a
    circular block bootstrap that destroys any trend while preserving
    short-range dependence. A slope outside the band is distinguishable
    from zero."""
    v = np.asarray(s, float)
    x = np.asarray(s.index, float)
    xc = x - x.mean()
    den = (xc ** 2).sum()
    n = len(v)
    slope = (v @ xc) / den
    starts = rng.integers(0, n, size=(nb, int(np.ceil(n / L))))
    idx = (starts[:, :, None] + np.arange(L)[None, None, :]).reshape(nb, -1)[:, :n] % n
    lo, hi = np.percentile((v[idx] @ xc) / den, [2.5, 97.5])
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
