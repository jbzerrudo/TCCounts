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
            ok("PAR mean per year, 1991-2020 standard normal",
               tot.loc[1991:2020].mean(), 20.2, .05)
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

print("\n--- Section 2, overlap between requirements ---")
for dom, sl, vd in [("WNP basin", 73, 199), ("PAR", 66, 197)]:
    a = R[(R.domain == dom) & (R.agency == "raw total")].sort_values(["rule", "start"])
    b = R[(R.domain == dom) & (R.agency == "all four")].sort_values(["rule", "start"])
    ok(f"{dom}: raw total and all four, identical slopes",
       (np.abs(a.slope.values - b.slope.values) < 1e-12).sum(), sl)
    ok(f"{dom}: raw total and all four, identical verdicts",
       (a.verdict.values == b.verdict.values).sum(), vd)

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

print("\n--- Section 2, dispersion and serial dependence of the counts ---")


def _fit(c):
    """The estimator from multiverse.py: slope and verdict at the 5% level."""
    v = np.asarray(c, float)
    x = np.asarray(c.index, float)
    xc = x - x.mean()
    den = (xc ** 2).sum()
    n = len(v)
    slope = (v @ xc) / den
    resid = v - np.poly1d(np.polyfit(x, v, 1))(x)
    r1 = min(max(np.corrcoef(resid[:-1], resid[1:])[0, 1], 0.0), 0.99)
    ne = max(n * (1 - r1) / (1 + r1), 4.0)
    se = np.sqrt((resid @ resid) / (n - 2) / den) * np.sqrt(n / ne)
    from scipy import stats as _st
    crit = _st.t.ppf(0.975, max(ne - 2, 1))
    return slope, ("+" if slope > crit * se else "-" if slope < -crit * se else "0")


def _r1(c):
    v = np.asarray(c, float)
    x = np.asarray(c.index, float)
    r = v - np.poly1d(np.polyfit(x, v, 1))(x)
    return np.corrcoef(r[:-1], r[1:])[0, 1]


_DOM = {}
if PAR:
    _DOM["PAR"] = pd.read_csv(PAR, low_memory=False)
if WNP:
    _DOM["WNP basin"] = pd.read_csv(WNP, low_memory=False, skiprows=[1])

for _dom, _raw in _DOM.items():
    _raw = _raw[(_raw.SEASON >= 1884) & (_raw.SEASON <= 2023)]
    _mn = _raw[_raw.TRACK_TYPE == "main"]
    _gc = "TOKYO_GRADE" if "TOKYO_GRADE" in _mn.columns else "TOK_GRADE"
    _msk = {"JMA": num(_mn, _gc).isin([2, 3, 4, 5, 9]), "JTWC": num(_mn, "USA_WIND") > 0,
            "CMA": num(_mn, "CMA_WIND") > 0, "HKO": num(_mn, "HKO_WIND") > 0}
    _tot = counts(_mn).loc[1951:2023]
    _v = np.asarray(_tot, float)
    ok(f"{_dom}: variance-to-mean ratio of the counts 1951-2023",
       _v.var(ddof=1) / _v.mean(), 1.37 if _dom == "WNP basin" else 1.01, .006)
    ok(f"{_dom}: lag-one autocorrelation of the detrended counts",
       _r1(_tot), 0.33 if _dom == "WNP basin" else 0.03, .006)

print("\n--- Section 3.4 and 4, single-agency counts and the coverage mechanism ---")
_WANT_SLOPE = {"WNP basin": {"CMA": (-0.166, "-"), "HKO": (0.203, "0"), "JTWC": (0.077, "0")},
               "PAR": {"CMA": (-0.091, "-"), "HKO": (0.135, "0"), "JTWC": (0.044, "+")}}
_WANT_COV = {"WNP basin": {"JTWC": 0.59, "HKO": 0.97, "CMA": -0.04},
             "PAR": {"JTWC": 0.59, "HKO": 0.94, "CMA": 0.02}}
_WANT_DECOMP = {"WNP basin": {"JTWC": (-0.137, 0.208), "HKO": (-0.120, 0.338), "CMA": (-0.152, -0.015)}}
_WANT_NSTART = {"WNP basin": {"CMA": 37}, "PAR": {"CMA": 26}}
_WANT_INTER = {"WNP basin": (-0.010, 0.196), "PAR": (-0.004, 0.133)}

for _dom, _raw in _DOM.items():
    _raw = _raw[(_raw.SEASON >= 1884) & (_raw.SEASON <= 2023)]
    _mn = _raw[_raw.TRACK_TYPE == "main"]
    _gc = "TOKYO_GRADE" if "TOKYO_GRADE" in _mn.columns else "TOK_GRADE"
    _msk = {"JMA": num(_mn, _gc).isin([2, 3, 4, 5, 9]), "JTWC": num(_mn, "USA_WIND") > 0,
            "CMA": num(_mn, "CMA_WIND") > 0, "HKO": num(_mn, "HKO_WIND") > 0}
    _totall = counts(_mn)
    _T = np.asarray(_totall.loc[1951:2023], float)
    _x = np.arange(1951, 2024, dtype=float)
    _bT = np.polyfit(_x, _T, 1)[0]
    for _a in ("JTWC", "CMA", "HKO"):
        _c = counts(_mn, _msk[_a])
        _sl, _vd = _fit(_c.loc[1951:2023])
        if _a in _WANT_SLOPE[_dom]:
            _w, _wv = _WANT_SLOPE[_dom][_a]
            ok(f"{_dom}, {_a} alone, slope 1951-2023", _sl, _w, .0006)
            ok(f"{_dom}, {_a} alone, verdict 1951-2023", _vd == _wv, True)
        _sh = 100 * np.asarray(_c.loc[1951:2023], float) / _T
        ok(f"{_dom}: {_a} coverage trend, percentage points per year",
           np.polyfit(_x, _sh, 1)[0], _WANT_COV[_dom][_a], .006)
        if _dom in _WANT_DECOMP:
            _wt, _wc = _WANT_DECOMP[_dom][_a]
            _s = np.asarray(_c.loc[1951:2023], float) / _T
            ok(f"{_dom}: {_a} count trend, term from the falling total",
               _s.mean() * _bT, _wt, .0006)
            ok(f"{_dom}: {_a} count trend, term from the rising share",
               _T.mean() * np.polyfit(_x, _s, 1)[0], _wc, .0006)
        if _a in _WANT_NSTART[_dom]:
            _n = sum(_fit(_c.loc[y:2023])[1] == "-" for y in range(1951, 2001))
            ok(f"{_dom}: {_a} alone, start years with a significant decline",
               _n, _WANT_NSTART[_dom][_a])
    _sid = {k: set(_mn[v].SID.unique()) for k, v in _msk.items()}
    _i3 = _sid["JMA"] & _sid["JTWC"] & _sid["CMA"]
    _i4 = _i3 & _sid["HKO"]
    for _lab, _ss, _w in [("JMA, JTWC and CMA", _i3, _WANT_INTER[_dom][0]),
                          ("all four including HKO", _i4, _WANT_INTER[_dom][1])]:
        ok(f"{_dom} intersection, {_lab}, slope 1951-2023",
           _fit(counts(_mn, _mn.SID.isin(_ss)).loc[1951:2023])[0], _w, .0006)
    ok(f"{_dom} four-way intersection, zero years 1951-1960",
       (counts(_mn, _mn.SID.isin(_i4)).loc[1951:1960] == 0).sum(), 10)

print("\n--- Section 4, the longest window and what survives the strictest exclusion ---")
for _dom, _dec, _non, _inc in [("WNP basin", 12, 12, 0), ("PAR", 12, 11, 1)]:
    _g = R[(R.domain == _dom) & (R.start == 1951)]
    ok(f"{_dom}: specifications beginning in 1951", len(_g), 24)
    ok(f"{_dom}: of those, significant declines", (_g.verdict == "-").sum(), _dec)
    ok(f"{_dom}: of those, no detectable trend", (_g.verdict == "0").sum(), _non)
    ok(f"{_dom}: of those, significant increases", (_g.verdict == "+").sum(), _inc)
    _cma = _g[_g.agency.isin(["JMA+JTWC+CMA", "all four", "raw total"])]
    _jj = _g[_g.agency.isin(["JMA", "JMA+JTWC", "JTWC"])]
    ok(f"{_dom}: at 1951, requirements admitting CMA or no grade that decline",
       (_cma.verdict == "-").sum(), 12)
    ok(f"{_dom}: at 1951, requirements restricted to JMA and JTWC that decline",
       (_jj.verdict == "-").sum(), 0)
    _rt = R[(R.domain == _dom) & (R.agency == "raw total") & (R.rule == "main track") & (R.start == 1951)].slope.iloc[0]
    _af = R[(R.domain == _dom) & (R.agency == "all four") & (R.rule == "main track") & (R.start == 1951)].slope.iloc[0]
    ok(f"{_dom}: share of the slope surviving the removal of ungraded storms (%)",
       100 * _af / _rt, 80 if _dom == "WNP basin" else 71, 0.6)
ok("windows beginning after 1990 are shorter than 34 years", 2023 - 1991 + 1, 33)

print("\n--- Section 2, the study-area map ---")
if WNP:
    from matplotlib.path import Path as _MPath   # same test clip_par.py uses
    _V = [(115, 5), (115, 15), (120, 21), (120, 25), (135, 25), (135, 5)]
    _hex = _MPath(_V)
    _b = pd.read_csv(WNP, low_memory=False, skiprows=[1],
                     usecols=["SID", "SEASON", "LAT", "LON", "TRACK_TYPE"])
    _b = _b[(_b.SEASON >= 1884) & (_b.SEASON <= 2023) & (_b.TRACK_TYPE == "main")].copy()
    _b["LAT"] = pd.to_numeric(_b.LAT, errors="coerce")
    _b["LON"] = pd.to_numeric(_b.LON, errors="coerce")
    _b = _b.dropna(subset=["LAT", "LON"])
    _box = _b[(_b.LAT >= 5) & (_b.LAT <= 25) & (_b.LON >= 115) & (_b.LON <= 135)]
    _in = _hex.contains_points(np.c_[_box.LON.values, _box.LAT.values], radius=-1e-9)
    _bs, _hs = set(_box.SID.unique()), set(_box.SID.values[_in])
    ok("storms a bounding box admits that the PAR excludes", len(_bs - _hs), 129)
    ok("storms inside the PAR hexagon, cross-check", len(_hs), 2757)
    ok("that overshoot as a percentage of the PAR count",
       100 * len(_bs - _hs) / len(_hs), 4.7, .05)
    _lon = _b.LON.to_numpy(float) % 360
    _lat = _b.LAT.to_numpy(float)
    _fx, _fy = np.floor(_lon - 100).astype(int), np.floor(_lat - 0).astype(int)
    _fr = (_fx >= 0) & (_fx < 80) & (_fy >= 0) & (_fy < 55)   # same binning as the figure
    ok("percentage of basin fixes inside the map frame", 100 * _fr.mean(), 97.5, .05)
    _pr = pd.DataFrame({"c": _fx[_fr] * 1000 + _fy[_fr],
                        "s": _b.SID.to_numpy()[_fr]}).drop_duplicates()
    ok("storms in the densest 1-degree cell of the map", _pr.groupby("c").size().max(), 196)

print("\n%d checks run, %d failed." % (N[0], FAIL[0]))
sys.exit(1 if FAIL[0] else 0)
