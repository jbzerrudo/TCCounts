"""
Version 4, step 1. Per-storm flags and the threshold grid.

Adds an intensity axis to the specification grid (all systems, tropical
depression or stronger, tropical storm or stronger, each agency judged by its
own categories; see classes.py) and rebuilds the grading requirements around
the three basin-wide agencies. HKO alone and any of the four agencies are kept
as separate series. Every series is tested with multiverse.fit, unchanged.

Outputs
    data/storm_flags.csv     one row per storm, domain and counting rule
    data/grid_threshold.csv  one row per specification (10,000 rows)
    data/archive_meta.csv    TRACK_TYPE values in the archive, by period
    data/cma_consistency.csv how often CMA_CAT disagrees with CMA_WIND

At the 'all' threshold, the requirements JMA, JTWC, JMA|JTWC, any of 3 and
any of 4, and the raw total, reproduce data/multiverse.csv exactly
(JMA+JTWC, JMA+JTWC+CMA and all four there); check_v4.py tests this.

Usage
    python src/threshold.py ibtracs.WP.list.v04r01.csv par_clipped.csv data/
"""
import os
import sys
import warnings

import pandas as pd

from classes import FLAGS, REQUIREMENTS, THRESHOLDS, annual, any4, requirement, storm_flags
from multiverse import FIRST, LAST, STARTS, counting_rules, fit


warnings.filterwarnings("ignore", category=RuntimeWarning)  # np.corrcoef on all-zero series


def main():
    wnp = sys.argv[1] if len(sys.argv) > 1 else "ibtracs.WP.list.v04r01.csv"
    par = sys.argv[2] if len(sys.argv) > 2 else "par_clipped.csv"
    out = sys.argv[3] if len(sys.argv) > 3 else "data"

    basin = pd.read_csv(wnp, low_memory=False, skiprows=[1])
    clip = pd.read_csv(par, low_memory=False)

    tt = basin.TRACK_TYPE.astype(str)
    meta = pd.DataFrame({"n_1884_2023": tt[(basin.SEASON >= FIRST) & (basin.SEASON <= LAST)].value_counts(),
                         "n_after_2023": tt[basin.SEASON > LAST].value_counts()}).fillna(0).astype(int)
    meta.index.name = "TRACK_TYPE"
    meta.to_csv(os.path.join(out, "archive_meta.csv"))

    # How often CMA's category disagrees with its own wind (fixes, 1884-2023)
    b = basin[(basin.SEASON >= FIRST) & (basin.SEASON <= LAST)]
    cw = pd.to_numeric(b.CMA_WIND, errors="coerce")
    cc = pd.to_numeric(b.CMA_CAT, errors="coerce")
    pd.DataFrame([dict(fixes_with_cma_wind=int((cw > 0).sum()),
                       ts_category_wind_below_33kt=int((cc.between(2, 6) & (cw < 33)).sum()),
                       td_or_weaker_category_wind_34kt_plus=int((cc.isin([0, 1]) & (cw >= 34)).sum()))]
                 ).to_csv(os.path.join(out, "cma_consistency.csv"), index=False)

    flag_rows, rows = [], []
    for dom, d in [("WNP basin", basin), ("PAR", clip)]:
        d = d[(d.SEASON >= FIRST) & (d.SEASON <= LAST)].copy()
        for rule, sub in counting_rules(d).items():
            fl = storm_flags(sub)
            f = fl.reset_index()
            f.insert(0, "rule", rule)
            f.insert(0, "domain", dom)
            flag_rows.append(f)
            series = {("raw total", "-"): annual(fl, pd.Series(True, index=fl.index))}
            for t in THRESHOLDS:
                for req in REQUIREMENTS:
                    series[(req, t)] = annual(fl, requirement(fl, req, t))
                series[("HKO", t)] = annual(fl, fl[f"HKO_{t}"])
                series[("any of 4", t)] = annual(fl, any4(fl, t))
            for (req, t), s in series.items():
                for y in STARTS:
                    slope, lo, hi, verdict = fit(s.loc[y:LAST])
                    rows.append(dict(domain=dom, rule=rule, requirement=req, threshold=t, start=y,
                                     slope=slope, lo=lo, hi=hi, verdict=verdict))
    flags = pd.concat(flag_rows, ignore_index=True)
    for c in FLAGS:
        flags[c] = flags[c].astype(int)
    flags.to_csv(os.path.join(out, "storm_flags.csv"), index=False)
    grid = pd.DataFrame(rows)
    grid.to_csv(os.path.join(out, "grid_threshold.csv"), index=False)

    for dom, g in grid.groupby("domain", sort=False):
        for t in THRESHOLDS:
            h = g[(g.threshold == t) & g.requirement.isin(REQUIREMENTS)]
            print(f"{dom:10s} {t:3s} n={len(h)}  significant declines {(h.verdict == '-').mean():6.1%}  "
                  f"median slope {h.slope.median():+.3f}")
    print(f"written: {out}/storm_flags.csv ({len(flags)} rows), {out}/grid_threshold.csv ({len(grid)} rows), "
          f"{out}/archive_meta.csv")


if __name__ == "__main__":
    main()
