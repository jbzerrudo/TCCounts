"""
Step 2. Build the three annual PAR count series from the clipped archive.

    total      every unique storm with a fix inside the PAR
    any_int    CLASSIFIED: some agency reported an intensity at some PAR fix,
               meaning a JMA grade of Tropical Depression or stronger, or a JTWC wind
    pos_only   UNCLASSIFIED: neither, so the storm is a position on a map and nothing else
    jma_class  the stricter JMA-graded-only subset, reported for sensitivity

total = any_int + pos_only exactly, every year. The two classes partition the total,
so their fitted trends sum to the fitted trend of the total. That identity is the
whole argument: it makes the decomposition an accounting statement, not a model.

Usage   python src/build_series.py par_clipped.csv data/par_annual_series.csv
"""
import sys
import pandas as pd

FIRST, LAST = 1884, 2023


def build(df, first=FIRST, last=LAST):
    # Spur tracks are excluded. IBTrACS gives a diverging agency track its own SID,
    # so counting every SID double-counts those storms; the IBTrACS documentation
    # lists spurs among the reasons to exercise care when counting. Forty-three PAR
    # identifiers are spur-only.
    d = df[(df.SEASON >= first) & (df.SEASON <= last) & (df.TRACK_TYPE == "main")].copy()
    gcol = "TOKYO_GRADE" if "TOKYO_GRADE" in d.columns else "TOK_GRADE"
    grade = pd.to_numeric(d[gcol], errors="coerce").fillna(0)
    wind = pd.to_numeric(d["USA_WIND"], errors="coerce").fillna(0)
    # JMA grades: 2 TD, 3 TS, 4 STS, 5 TY, 9 TC of TS intensity or higher.
    # 6 is extratropical and 7 is a position marker, so neither is an intensity
    # classification under the definition used in the paper.
    d["_jma"] = grade.isin([2, 3, 4, 5, 9])
    d["_any"] = d["_jma"] | (wind > 0)
    years = range(first, last + 1)
    n = lambda sub: sub.groupby("SEASON").SID.nunique().reindex(years).fillna(0).astype(int)
    s = pd.DataFrame({"total": n(d), "jma_class": n(d[d._jma]), "any_int": n(d[d._any])})
    s["pos_only"] = s.total - s.any_int
    s.index.name = "SEASON"
    return s


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "par_clipped.csv"
    dst = sys.argv[2] if len(sys.argv) > 2 else "data/par_annual_series.csv"
    s = build(pd.read_csv(src, low_memory=False))
    s.to_csv(dst)
    print(s.describe().round(2).to_string())
    print(f"written: {dst}")


if __name__ == "__main__":
    main()
