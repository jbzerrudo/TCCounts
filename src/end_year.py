"""
Version 4, step 6. Does the verdict depend on the end year?

Main track. For end years 2005, 2010, 2015, 2020 and 2023, every window that
starts in 1951 or later and is at least 24 years long is fitted with the test
in multiverse.py; the share of windows with a significant decline and the
slope from 1951 are recorded.

Output
    data/end_year.csv

Usage
    python src/end_year.py data/storm_flags.csv data/
"""
import os
import sys

import numpy as np
import pandas as pd

from classes import annual, fit_many, load_flags, requirement

ENDS = (2005, 2010, 2015, 2020, 2023)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "data/storm_flags.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "data"
    rows = []
    for dom in ("WNP basin", "PAR"):
        fl = load_flags(src, dom)
        ser = {"raw total": pd.Series(True, index=fl.index),
               "any of 3, all systems": requirement(fl, "any of 3", "all"),
               "JMA|JTWC, all systems": requirement(fl, "JMA|JTWC", "all"),
               "JMA, TS": fl.JMA_TS, "JTWC, TS": fl.JTWC_TS, "CMA, TS": fl.CMA_TS,
               "any of 3, TS": requirement(fl, "any of 3", "TS"),
               "all 3, TS": requirement(fl, "all 3", "TS")}
        for name, sel in ser.items():
            s = annual(fl, sel)
            for e in ENDS:
                starts = np.arange(1951, e - 23 + 1)
                dec, first = 0, None
                for y in starts:
                    sl, band, vd = fit_many(s.loc[y:e].to_numpy(float)[None, :], np.arange(y, e + 1))
                    dec += int(vd[0] == -1)
                    if y == 1951:
                        first = (sl[0], {-1: "-", 0: "0", 1: "+"}[int(vd[0])])
                rows.append(dict(domain=dom, series=name, end=e, windows=len(starts),
                                 share_sig_decline_pct=100 * dec / len(starts),
                                 slope_from_1951=first[0], verdict_from_1951=first[1]))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(out, "end_year.csv"), index=False)
    t = R.assign(cell=R.share_sig_decline_pct.map(lambda v: f"{v:3.0f}%") + " | " +
                 R.slope_from_1951.map(lambda v: f"{v:+.3f}") + R.verdict_from_1951.map({"-": "*", "+": "^", "0": " "}))
    pd.set_option("display.width", 200)
    print("share of windows with a significant decline | slope from 1951 (* decline, ^ increase)")
    print(t.pivot_table(index=["domain", "series"], columns="end", values="cell", aggfunc="first").to_string())
    print(f"written: {out}/end_year.csv")


if __name__ == "__main__":
    main()
