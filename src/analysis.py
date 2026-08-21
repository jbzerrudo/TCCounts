"""
Recompute every value quoted in the manuscript and check each one.
Prints OK or MISMATCH beside every value. A clean run prints zero MISMATCH.

Usage
    python src/analysis.py data/multiverse.csv par_clipped.csv \
                           ibtracs.WP.list.v04r01.csv
"""
import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

MV = sys.argv[1] if len(sys.argv) > 1 else "data/multiverse.csv"
PAR = sys.argv[2] if len(sys.argv) > 2 else "par_clipped.csv"
WNP = sys.argv[3] if len(sys.argv) > 3 else None

R = pd.read_csv(MV)
N, FAIL = [0], [0]
AGS = ["raw total", "all four", "JMA+JTWC+CMA", "JMA", "JMA+JTWC", "JTWC"]
RUS = ["all SIDs", "main track", ">=2 fixes", "tropical nature"]
YEARS = list(range(1884, 2024))


def ok(name, got, want, tol=0):
    N[0] += 1
    good = abs(float(got) - want) <= tol
    if not good:
        FAIL[0] += 1
    print(("  OK   " if good else "  MISMATCH") +
          "  %-48s paper=%-10s computed=%s" % (name, want, round(float(got), 4)))


def num(d, c):
    return pd.to_numeric(d[c], errors="coerce").fillna(0)


def counts(d, mask=None):
    sub = d if mask is None else d[mask]
    return sub.groupby("SEASON").SID.nunique().reindex(YEARS).fillna(0).astype(int)


def pct_decline(g):
    return 100 * (g.verdict == "-").mean()


print("\n--- Section 2, the specification space ---")
ok("specifications per domain", len(R[R.domain == "PAR"]), 1200)
ok("specifications in total", len(R), 2400)
ok("agency sets", R.agency.nunique(), 6)
ok("counting rules", R.rule.nunique(), 4)
ok("start years", R.start.nunique(), 50)
ok("first start year", R.start.min(), 1951)
ok("last start year", R.start.max(), 2000)

print("\n--- Section 3.1, how the space divides ---")
ok("significant increases, all domains", (R.verdict == "+").sum(), 0)
for dom, dec, non, lo, hi in [("WNP basin", 498, 702, -0.3488, 0.0964),
                              ("PAR", 364, 836, -0.1595, 0.1062)]:
    g = R[R.domain == dom]
    ok(f"{dom}: significant declines", (g.verdict == "-").sum(), dec)
    ok(f"{dom}: no detectable trend", (g.verdict == "0").sum(), non)
    ok(f"{dom}: percent declining", pct_decline(g), round(100 * dec / 1200), 0.5)
    ok(f"{dom}: minimum slope", g.slope.min(), lo, .0005)
    ok(f"{dom}: maximum slope", g.slope.max(), hi, .0005)

print("\n--- Section 3.1, the two worked specifications ---")
for dom, a, want, verd in [("PAR", "raw total", -0.1086, "-"), ("PAR", "JMA+JTWC", 0.0131, "0"),
                           ("WNP basin", "raw total", -0.1644, "-"),
                           ("WNP basin", "JMA+JTWC", 0.0439, "0")]:
    r = R[(R.domain == dom) & (R.agency == a) & (R.rule == "main track") & (R.start == 1951)].iloc[0]
    ok(f"{dom}, {a}, main track, from 1951: slope", r.slope, want, .0005)
    ok(f"{dom}, {a}, main track, from 1951: verdict", r.verdict == verd, True)

print("\n--- Section 3.2, agency choice dominates ---")
for dom, vals in [("WNP basin", [68, 68, 67, 18, 12, 14]), ("PAR", [59, 54, 56, 10, 0, 2])]:
    for a, v in zip(AGS, vals):
        ok(f"{dom}: {a}, percent declining", pct_decline(R[(R.domain == dom) & (R.agency == a)]), v, 0.5)
for dom, a, v in [("WNP basin", "raw total", -0.197), ("WNP basin", "JMA+JTWC", -0.019),
                  ("PAR", "raw total", -0.114), ("PAR", "JMA+JTWC", -0.006)]:
    ok(f"{dom}: {a}, median slope", R[(R.domain == dom) & (R.agency == a)].slope.median(), v, .0006)

print("\n--- Section 3.2, the counting rule is not a third knob ---")
for dom, vals in [("WNP basin", [40, 42, 41, 43]), ("PAR", [36, 30, 32, 24])]:
    for r_, v in zip(RUS, vals):
        ok(f"{dom}: {r_}, percent declining", pct_decline(R[(R.domain == dom) & (R.rule == r_)]), v, 0.5)

print("\n--- Section 3.3, the start year is the second knob ---")
R["decade"] = (R.start // 10) * 10
for dom, vals in [("PAR", [52, 57, 38, 9, 2, 0]), ("WNP basin", [52, 62, 18, 31, 49, 0])]:
    for (dec, g), v in zip(R[R.domain == dom].groupby("decade"), vals):
        ok(f"{dom}: start years in the {int(dec)}s, percent declining", pct_decline(g), v, 0.5)

if PAR:
    print("\n--- Section 3.4, agency coverage ---")
    doms = [("PAR", pd.read_csv(PAR, low_memory=False))]
    if WNP:
        doms.append(("WNP basin", pd.read_csv(WNP, low_memory=False, skiprows=[1])))
    want = {"PAR": {"JMA": [71, 77, 84], "JTWC": [67, 88, 96], "CMA": [90, 91, 93], "HKO": [45, 78, 85]},
            "WNP basin": {"JMA": [68, 77, 80], "JTWC": [67, 89, 96], "CMA": [93, 91, 92], "HKO": [50, 85, 89]}}
    eras = [(1951, 1976), (1977, 2000), (2001, 2023)]
    for dom, d in doms:
        d = d[(d.SEASON >= 1884) & (d.SEASON <= 2023)]
        m = d[d.TRACK_TYPE == "main"]
        gcol = "TOKYO_GRADE" if "TOKYO_GRADE" in m.columns else "TOK_GRADE"
        tot = counts(m)
        masks = {"JMA": num(m, gcol).isin([2, 3, 4, 5, 9]), "JTWC": num(m, "USA_WIND") > 0,
                 "CMA": num(m, "CMA_WIND") > 0, "HKO": num(m, "HKO_WIND") > 0}
        for a, vals in want[dom].items():
            c = counts(m, masks[a])
            for (x, y), v in zip(eras, vals):
                ok(f"{dom}: {a} coverage {x}-{y} (%)",
                   100 * c.loc[x:y].sum() / tot.loc[x:y].sum(), v, 0.6)
        if dom == "PAR":
            ok("PAR storms 1884-2023", tot.sum(), 2757)
            ok("PAR mean per year", tot.mean(), 19.7, .05)
            ok("PAR storms 1884-1944", tot.loc[1884:1944].sum(), 1004)
            any_int = counts(m, masks["JMA"] | masks["JTWC"] | masks["CMA"] | masks["HKO"])
            ok("PAR storms 1884-1944 with any agency intensity", any_int.loc[1884:1944].sum(), 0)
            for a, y in [("JTWC", 1945), ("CMA", 1949), ("JMA", 1951), ("HKO", 1961)]:
                first = counts(m, masks[a])
                ok(f"PAR: first season {a} reports an intensity",
                   first[first > 0].index.min(), y)
        else:
            ok("WNP basin storms 1884-2023", tot.sum(), 4076)
            ok("WNP basin mean per year", tot.mean(), 29.1, .05)

print("\n%d checks run, %d failed." % (N[0], FAIL[0]))
sys.exit(1 if FAIL[0] else 0)
