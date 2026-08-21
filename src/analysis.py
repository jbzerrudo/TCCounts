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
ok("significant increases, all domains", (R.verdict == "+").sum(), 1)
for dom, dec, non, inc, pct, lo, hi in [("WNP basin", 620, 580, 0, 51.7, -0.3488, 0.0964),
                                        ("PAR", 400, 799, 1, 33.3, -0.1595, 0.1062)]:
    g = R[R.domain == dom]
    ok(f"{dom}: significant declines", (g.verdict == "-").sum(), dec)
    ok(f"{dom}: no detectable trend", (g.verdict == "0").sum(), non)
    ok(f"{dom}: significant increases", (g.verdict == "+").sum(), inc)
    ok(f"{dom}: percent declining", pct_decline(g), pct, .05)
    ok(f"{dom}: minimum slope", g.slope.min(), lo, .0005)
    ok(f"{dom}: maximum slope", g.slope.max(), hi, .0005)
    sub = g[g.agency != "raw total"]
    ok(f"{dom}: percent declining, raw total excluded", pct_decline(sub),
       44.8 if dom == "WNP basin" else 27.3, .05)

inc = R[R.verdict == "+"].iloc[0]
ok("the increase is in the PAR", inc.domain == "PAR", True)
ok("the increase restricts the count to JTWC", inc.agency == "JTWC", True)
ok("the increase uses the main-track rule", inc.rule == "main track", True)
ok("the increase begins in 1951", inc.start, 1951)
ok("increase slope, PAR main track JTWC 1951", inc.slope, 0.044428, .0005)
ok("PAR, JMA or JTWC, significant declines",
   ((R.domain == "PAR") & (R.agency == "JMA+JTWC") & (R.verdict == "-")).sum(), 0)

print("\n--- Section 3.1, the two worked specifications ---")
for dom, a, want, verd in [("PAR", "raw total", -0.1086, "-"), ("PAR", "JMA+JTWC", 0.0131, "0"),
                           ("WNP basin", "raw total", -0.1644, "-"),
                           ("WNP basin", "JMA+JTWC", 0.0439, "0")]:
    r = R[(R.domain == dom) & (R.agency == a) & (R.rule == "main track") & (R.start == 1951)].iloc[0]
    ok(f"{dom}, {a}, main track, from 1951: slope", r.slope, want, .0005)
    ok(f"{dom}, {a}, main track, from 1951: verdict", r.verdict == verd, True)

print("\n--- Section 3.2, agency choice dominates ---")
for dom, vals in [("WNP basin", [86.0, 85.5, 87.0, 27.5, 12.5, 11.5]),
                  ("PAR", [63.5, 62.0, 62.5, 12.0, 0.0, 0.0])]:
    for a, v in zip(AGS, vals):
        ok(f"{dom}: {a}, percent declining", pct_decline(R[(R.domain == dom) & (R.agency == a)]), v, .05)
for dom, vals in [("WNP basin", [-0.197, -0.176, -0.174, -0.086, -0.019, -0.005]),
                  ("PAR", [-0.114, -0.103, -0.101, -0.041, -0.014, -0.007])]:
    for a, v in zip(AGS, vals):
        ok(f"{dom}: {a}, median slope",
           R[(R.domain == dom) & (R.agency == a)].slope.median(), v, .0006)

print("\n--- Section 3.2, the counting rule is not a third knob ---")
for dom, vals in [("WNP basin", [56.3, 53.0, 52.7, 44.7]), ("PAR", [37.0, 30.3, 34.7, 31.3])]:
    for r_, v in zip(RUS, vals):
        ok(f"{dom}: {r_}, percent declining", pct_decline(R[(R.domain == dom) & (R.rule == r_)]), v, .05)

print("\n--- Section 3.3, the start year is the second knob ---")
R["decade"] = (R.start // 10) * 10
for dom, vals in [("PAR", [56, 55, 49, 12, 0, 0]), ("WNP basin", [60, 63, 49, 57, 36, 0])]:
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

print("\n--- Section 2 and 3.2, structure of the specification space ---")
for dom, vals in [("WNP basin", (40.7, 29.1, 2.9)), ("PAR", (38.3, 33.3, 0.7))]:
    g = R[R.domain == dom]
    tot = g.slope.var(ddof=0)
    for f, want in zip(["agency", "start", "rule"], vals):
        share = ((g.groupby(f).slope.mean() - g.slope.mean()) ** 2 * g.groupby(f).size()).sum() / len(g)
        ok(f"{dom}: {f} share of slope variance (%)", 100 * share / tot, want, .05)
    M = g[g.rule == "main track"].pivot(index="start", columns="agency", values="slope").corr().values
    ev = np.linalg.eigvalsh(M)[::-1]
    ok(f"{dom}: effective dimensionality of the six requirements",
       (ev.sum() ** 2) / (ev ** 2).sum(), 1.8, .05)

print("\n--- Section 4, the intersection sets ---")
if PAR:
    _d = pd.read_csv(PAR, low_memory=False)
    _d = _d[(_d.SEASON >= 1884) & (_d.SEASON <= 2023)]
    _m = _d[_d.TRACK_TYPE == "main"]
    _g = "TOKYO_GRADE" if "TOKYO_GRADE" in _m.columns else "TOK_GRADE"
    _per = lambda k: set(_m[k].SID.unique())
    _J = _per(num(_m, _g).isin([2, 3, 4, 5, 9])); _T = _per(num(_m, "USA_WIND") > 0)
    _C = _per(num(_m, "CMA_WIND") > 0); _H = _per(num(_m, "HKO_WIND") > 0)
    for lab, ss, want in [("JMA, JTWC and CMA", _J & _T & _C, -0.004),
                          ("all four including HKO", _J & _T & _C & _H, 0.133)]:
        c = counts(_m, _m.SID.isin(ss)).loc[1951:2023]
        x = np.asarray(c.index, float); xc = x - x.mean()
        ok(f"PAR intersection, {lab}, slope 1951-2023",
           (np.asarray(c, float) @ xc) / (xc ** 2).sum(), want, .0006)
    _hko = counts(_m, _m.SID.isin(_J & _T & _C & _H)).loc[1951:1960]
    ok("PAR four-way intersection, zero years 1951-1960", (_hko == 0).sum(), 10)

    print("\n--- Section 3.4, single-agency counts ---")
    _slope = lambda c: (np.asarray(c, float) @ (np.asarray(c.index, float) - np.asarray(c.index, float).mean())) \
                       / (((np.asarray(c.index, float) - np.asarray(c.index, float).mean())) ** 2).sum()
    for lab, msk, want in [("CMA alone", num(_m, "CMA_WIND") > 0, -0.091),
                           ("HKO alone", num(_m, "HKO_WIND") > 0, 0.135)]:
        ok(f"PAR, {lab}, slope 1951-2023", _slope(counts(_m, msk).loc[1951:2023]), want, .0006)
    _h = counts(_m, num(_m, "HKO_WIND") > 0)
    ok("PAR, first season HKO grades a storm", _h[_h > 0].index.min(), 1961)
    ok("PAR, HKO zero years 1951-1960", (_h.loc[1951:1960] == 0).sum(), 10)
    _k = num(_m, "KMA_WIND") > 0
    ok("PAR, first season KMA reports a wind", _m[_k].SEASON.min(), 2015)
    _int = [c for c in _m.columns if c.endswith(("_WIND", "_PRES", "_GRADE", "_CAT"))]
    ok("intensity-bearing columns in the archive", len(_int), 38)
    _pre = _m[(_m.SEASON >= 1884) & (_m.SEASON <= 1944)]
    ok("pre-1945 PAR fixes with any intensity, all 38 columns",
       sum(int((num(_pre, c) > 0).sum()) for c in _int), 0)
    ok("TRACK_TYPE values are main plus spurs only",
       set(_d.TRACK_TYPE.unique()) <= {"main", "spur-merge", "spur-split", "spur-other"}, True)

print("\n%d checks run, %d failed." % (N[0], FAIL[0]))
sys.exit(1 if FAIL[0] else 0)
