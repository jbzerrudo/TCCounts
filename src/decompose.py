"""
Version 4, step 2. Where the decline sits.

Each main-track storm is assigned to one class by the strongest category any of
the four agencies gives it inside the domain: tropical storm or stronger,
depression only, below depression (a positive wind but no depression category
from anyone), or no agency intensity. The four classes partition the total, and
OLS slopes add, so the split of the total trend is exact. The same is done by
which agencies graded the storm (all-systems definition).

Outputs
    data/class_counts.csv        annual counts per class, both domains, 1884-2023
    data/decomposition.csv       slope, verdict and share of the total slope per
                                 class, windows starting 1951, 1961 and 1977
    data/era_means.csv           storms per year per class, three eras
    data/ts_coverage.csv         share of tropical-storm-strength storms that each
                                 agency also rates tropical storm or stronger, by era
    data/cma_wind_sensitivity.csv  CMA classified by wind instead of category

Usage
    python src/decompose.py data/storm_flags.csv data/
"""
import os
import sys
import warnings

import pandas as pd

from classes import annual, any4, load_flags
from multiverse import LAST, fit

ERAS = [(1951, 1976), (1977, 2000), (2001, 2023)]
STARTS = (1951, 1961, 1977)


warnings.filterwarnings("ignore", category=RuntimeWarning)  # np.corrcoef on all-zero series


def classes(fl):
    ts, td, al = any4(fl, "TS"), any4(fl, "TD"), any4(fl, "all")
    ts3 = fl[["JMA_TS", "JTWC_TS", "CMA_TS"]]
    return {
        "total": pd.Series(True, index=fl.index),
        "ts_or_stronger": ts,
        "depression_only": td & ~ts,
        "below_depression": al & ~td,
        "no_intensity": ~al,
        "jma_graded": fl.JMA_all,
        "jtwc_not_jma": fl.JTWC_all & ~fl.JMA_all,
        "cma_hko_only": (fl.CMA_all | fl.HKO_all) & ~(fl.JMA_all | fl.JTWC_all),
        "hko_only": fl.HKO_all & ~(fl.JMA_all | fl.JTWC_all | fl.CMA_all),
        "ts_all3": ts3.all(axis=1),
        "ts_some3": ts3.any(axis=1) & ~ts3.all(axis=1),
        "ts_hko_only": fl.HKO_TS & ~ts3.any(axis=1),
    }


PARTITION = ("ts_or_stronger", "depression_only", "below_depression", "no_intensity")
BY_AGENCY = ("jma_graded", "jtwc_not_jma", "cma_hko_only", "no_intensity")


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "data/storm_flags.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "data"
    counts, dec, eras, cov, cma = [], [], [], [], []
    for dom in ("WNP basin", "PAR"):
        fl = load_flags(src, dom)
        cl = classes(fl)
        ser = {k: annual(fl, v) for k, v in cl.items()}
        c = pd.DataFrame(ser)
        c.index.name = "SEASON"
        c.insert(0, "domain", dom)
        counts.append(c.reset_index())
        for y in STARTS:
            tot = fit(ser["total"].loc[y:LAST])[0]
            for k in ser:
                s = ser[k].loc[y:LAST]
                sl, lo, hi, vd = fit(s)
                dec.append(dict(domain=dom, start=y, cls=k, slope=sl, verdict=vd,
                                share_of_total_pct=100 * sl / tot, mean_per_year=s.mean()))
        for k in ser:
            for a, b in ERAS:
                eras.append(dict(domain=dom, cls=k, era=f"{a}-{b}", mean_per_year=ser[k].loc[a:b].mean()))
        base = cl["ts_or_stronger"]
        for ag in ("JMA", "JTWC", "CMA", "HKO", "all three"):
            for a, b in ERAS:
                sel = base & fl.SEASON.between(a, b)
                if ag == "all three":
                    v = fl.loc[sel, ["JMA_TS", "JTWC_TS", "CMA_TS"]].all(axis=1).mean()
                else:
                    v = fl.loc[sel, f"{ag}_TS"].mean()
                cov.append(dict(domain=dom, agency=ag, era=f"{a}-{b}", pct=100 * v))
        J, T = fl.JMA_TS, fl.JTWC_TS
        sens = {"CMA, TS by category": fl.CMA_TS, "CMA, TS by wind": fl.CMA_TSwind,
                "any of 3, TS, CMA by category": J | T | fl.CMA_TS,
                "any of 3, TS, CMA by wind": J | T | fl.CMA_TSwind}
        for name, sel in sens.items():
            s = annual(fl, sel)
            for y in STARTS:
                sl, lo, hi, vd = fit(s.loc[y:LAST])
                cma.append(dict(domain=dom, series=name, start=y, slope=sl, verdict=vd))
    pd.concat(counts, ignore_index=True).to_csv(os.path.join(out, "class_counts.csv"), index=False)
    D = pd.DataFrame(dec)
    D.to_csv(os.path.join(out, "decomposition.csv"), index=False)
    pd.DataFrame(eras).to_csv(os.path.join(out, "era_means.csv"), index=False)
    pd.DataFrame(cov).to_csv(os.path.join(out, "ts_coverage.csv"), index=False)
    pd.DataFrame(cma).to_csv(os.path.join(out, "cma_wind_sensitivity.csv"), index=False)

    for dom in ("WNP basin", "PAR"):
        g = D[(D.domain == dom) & (D.start == 1951)].set_index("cls")
        print(f"{dom}, main track, 1951-2023: total {g.loc['total', 'slope']:+.4f} ({g.loc['total', 'verdict']})")
        for k in PARTITION:
            print(f"    {k:18s} {g.loc[k, 'slope']:+.4f} ({g.loc[k, 'verdict']})  "
                  f"{g.loc[k, 'share_of_total_pct']:6.1f}% of the total slope")
    print(f"written: class_counts, decomposition, era_means, ts_coverage, cma_wind_sensitivity in {out}/")


if __name__ == "__main__":
    main()
