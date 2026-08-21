"""
Build data/coverage.csv: the annual main-track storm count for each domain, and
how many of those storms each of the four agencies assigned an intensity to.

This is the series behind Fig 2c and 2d, Table 1, and the coverage decomposition
in Section 4. Counting follows multiverse.py exactly: main-track identifiers
only, one row per season and domain.

Usage
    python src/coverage.py ibtracs.WP.list.v04r01.csv par_clipped.csv \
                           data/coverage.csv
"""
import sys

import pandas as pd

from multiverse import FIRST, LAST, YEARS, num

AGENCIES = ("JMA", "JTWC", "CMA", "HKO")


def masks(m):
    g = "TOKYO_GRADE" if "TOKYO_GRADE" in m.columns else "TOK_GRADE"
    return {"JMA": num(m, g).isin([2, 3, 4, 5, 9]),
            "JTWC": num(m, "USA_WIND") > 0,
            "CMA": num(m, "CMA_WIND") > 0,
            "HKO": num(m, "HKO_WIND") > 0}


def main():
    wnp = sys.argv[1] if len(sys.argv) > 1 else "ibtracs.WP.list.v04r01.csv"
    par = sys.argv[2] if len(sys.argv) > 2 else "par_clipped.csv"
    dst = sys.argv[3] if len(sys.argv) > 3 else "data/coverage.csv"

    rows = []
    for domain, path in (("WNP basin", wnp), ("PAR", par)):
        d = pd.read_csv(path, low_memory=False,
                        skiprows=[1] if domain == "WNP basin" else None)
        d = d[(d.SEASON >= FIRST) & (d.SEASON <= LAST)]
        m = d[d.TRACK_TYPE == "main"]
        n = lambda sub: (sub.groupby("SEASON").SID.nunique()
                         .reindex(YEARS).fillna(0).astype(int))
        tot, mk = n(m), masks(m)
        cov = {a: n(m[mk[a]]) for a in AGENCIES}
        for y in YEARS:
            rows.append(dict(SEASON=y, domain=domain, total=int(tot.loc[y]),
                             **{a: int(cov[a].loc[y]) for a in AGENCIES}))
    out = pd.DataFrame(rows)
    out.to_csv(dst, index=False)
    print(f"{len(out)} rows written: {dst}")


if __name__ == "__main__":
    main()
