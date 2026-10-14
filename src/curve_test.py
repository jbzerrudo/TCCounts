"""
Version 4, step 5. A test of the specification curve as a whole.

Null hypothesis: no trend in any series. Synthetic archives are built by a
circular block bootstrap of whole years, 1951 to 2023; every series in a year
moves together, so the dependence between specifications is kept. Each
synthetic archive is run through every specification of a family, and the
number of significant declines and the median slope are recorded. The p-value
is the share of synthetic archives at least as extreme as the real one (in the
spirit of Simonsohn, Simmons and Nelson 2020, Nat Hum Behav 4:1208-1214).

Families, 24 series x 50 start years = 1,200 specifications per domain each:
    original  the six requirements of data/multiverse.csv (all-systems definition)
    all       the rebuilt requirements of classes.py, all systems
    TS        the rebuilt requirements at tropical-storm strength

Output
    data/curve_test.csv

Usage
    python src/curve_test.py data/storm_flags.csv data/ [replicates] [block lengths]
    defaults: 2000 replicates, block lengths 5,10,15 (about ten minutes)
"""
import os
import sys

import numpy as np
import pandas as pd

from classes import FLAGS, REQUIREMENTS, annual, fit_many, requirement

RULES = ("all SIDs", "main track", ">=2 fixes", "tropical nature")
STARTS = np.arange(1951, 2001)
YRS = np.arange(1951, 2024)


def families(flags, dom):
    fam = {"original": [], "all": [], "TS": []}
    for rule in RULES:
        fl = flags[(flags.domain == dom) & (flags.rule == rule)].set_index("SID")
        for c in FLAGS:
            fl[c] = fl[c].astype(bool)
        get = lambda sel: annual(fl, sel).loc[1951:2023].to_numpy(float)
        any4 = fl.JMA_all | fl.JTWC_all | fl.CMA_all | fl.HKO_all
        fam["original"] += [get(pd.Series(True, index=fl.index)), get(any4),
                            get(requirement(fl, "any of 3", "all")), get(fl.JMA_all),
                            get(requirement(fl, "JMA|JTWC", "all")), get(fl.JTWC_all)]
        fam["all"] += [get(requirement(fl, r, "all")) for r in REQUIREMENTS]
        fam["TS"] += [get(requirement(fl, r, "TS")) for r in REQUIREMENTS]
    return {k: np.vstack(v) for k, v in fam.items()}


def curve(M):
    dec = inc = 0
    slopes = []
    for y in STARTS:
        i = y - 1951
        sl, band, vd = fit_many(M[:, i:], YRS[i:])
        dec += int((vd == -1).sum())
        inc += int((vd == 1).sum())
        slopes.append(sl)
    return dec, inc, float(np.median(np.concatenate(slopes)))


def block_index(n, L, rng):
    st = rng.integers(0, n, size=int(np.ceil(n / L)))
    return ((st[:, None] + np.arange(L)[None, :]) % n).ravel()[:n]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "data/storm_flags.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "data"
    B = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
    LS = [int(v) for v in sys.argv[4].split(",")] if len(sys.argv) > 4 else [5, 10, 15]
    flags = pd.read_csv(src)
    rows = []
    for dom in ("WNP basin", "PAR"):
        F = families(flags, dom)
        obs = {k: curve(M) for k, M in F.items()}
        for L in LS:
            rng = np.random.default_rng(12345 + L)
            null = {k: [] for k in F}
            for _ in range(B):
                idx = block_index(73, L, rng)
                for k, M in F.items():
                    null[k].append(curve(M[:, idx]))
            for k in F:
                a = np.array(null[k])
                d0, i0, m0 = obs[k]
                rows.append(dict(domain=dom, family=k, block=L, replicates=B, specs=F[k].shape[0] * len(STARTS),
                                 obs_declines=d0, obs_increases=i0, obs_median_slope=m0,
                                 null_declines_median=float(np.median(a[:, 0])),
                                 null_declines_p95=float(np.percentile(a[:, 0], 95)),
                                 p_declines=float(np.mean(a[:, 0] >= d0)),
                                 p_median_slope=float(np.mean(a[:, 2] <= m0))))
            print(f"{dom:10s} block {L:2d}: done", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(out, "curve_test.csv"), index=False)
    pd.set_option("display.width", 200)
    print(R[["domain", "family", "block", "obs_declines", "null_declines_p95", "p_declines", "p_median_slope"]].to_string(index=False))
    print(f"written: {out}/curve_test.csv")


if __name__ == "__main__":
    main()
