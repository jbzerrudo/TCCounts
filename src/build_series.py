"""
Step 2. Build the three annual PAR count series from the clipped archive.

    total      every unique storm with a fix inside the PAR
    jma_class  storms the JMA graded Tropical Depression or stronger (grade 2-6)
    pos_only   storms with no JMA grade and no JTWC wind anywhere inside the PAR

total = classified-by-any-agency + pos_only, exactly, every year.

Usage   python src/build_series.py par_clipped.csv data/par_annual_series.csv
"""
import sys
import pandas as pd

FIRST, LAST = 1923, 2023


def build(df, first=FIRST, last=LAST):
    d = df[(df.SEASON >= first) & (df.SEASON <= last)].copy()
    grade = pd.to_numeric(d["TOK_GRADE"], errors="coerce").fillna(0)
    wind = pd.to_numeric(d["USA_WIND"], errors="coerce").fillna(0)
    d["_jma"] = grade.between(2, 6)
    d["_any"] = (grade > 0) | (wind > 0)
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
