"""
Recompute every value quoted in the manuscript and check each one.
Prints OK or MISMATCH beside every value. A clean run prints zero MISMATCH.

Usage   python src/analysis.py [par_annual_series.csv] [start_year_sweep.csv] [par_clipped.csv]

The third argument is optional and enables the counting-rule sensitivity checks
(Table 3), which need the clipped archive rather than the derived series.
"""
import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import warnings
warnings.filterwarnings("ignore")

SER = sys.argv[1] if len(sys.argv) > 1 else "data/par_annual_series.csv"
SWP = sys.argv[2] if len(sys.argv) > 2 else "data/start_year_sweep.csv"
CLIP = sys.argv[3] if len(sys.argv) > 3 else None

S = pd.read_csv(SER).set_index("SEASON")
sw = pd.read_csv(SWP)
rng = np.random.default_rng(41)
FAIL = [0]


N = [0]


def ok(name, got, want, tol):
    N[0] += 1
    good = abs(float(got) - want) <= tol
    if not good:
        FAIL[0] += 1
    print(("  OK   " if good else "  MISMATCH") + "  %-46s paper=%-10s computed=%s"
          % (name, want, round(float(got), 4)))


def slope(s):
    return st.linregress(np.asarray(s.index, float), np.asarray(s, float)).slope


def r2(s):
    return st.linregress(np.asarray(s.index, float), np.asarray(s, float)).rvalue ** 2


def band(s, L=10, nb=6000):
    v = np.asarray(s, float); x = np.asarray(s.index, float); N = len(v); out = []
    for _ in range(nb):
        i = np.concatenate([np.arange(k, k + L) % N
                            for k in rng.integers(0, N, size=int(np.ceil(N / L)))])[:N]
        out.append(st.linregress(x, v[i]).slope)
    return np.percentile(out, [2.5, 97.5])


print("\n--- Section 2, data and methods ---")
ok("annual values", len(S), 140, 0)
ok("storms counted", S.total.sum(), 2757, 0)
ok("mean per year", S.total.mean(), 19.69, .006)
ok("maximum", S.total.max(), 35, 0); ok("year of maximum", S.total.idxmax(), 1964, 0)
ok("minimum", S.total.min(), 8, 0); ok("year of minimum", S.total.idxmin(), 1885, 0)

print("\n--- Section 3.1, archive provenance ---")
ok("storms 1884-1944", S.loc[1884:1944, "total"].sum(), 1004, 0)
ok("classified 1884-1944", S.loc[1884:1944, "any_int"].sum(), 0, 0)
for a, b, v in [(1884, 1944, 100), (1945, 1950, 30), (1951, 1976, 26),
                (1977, 2000, 10), (2001, 2023, 2)]:
    ok("unclassified share %d-%d (%%)" % (a, b),
       100 * S.loc[a:b, "pos_only"].mean() / S.loc[a:b, "total"].mean(), v, .7)
for a, b, v in [(1884, 1922, 14.87), (1951, 1976, 25.12), (2001, 2023, 19.26)]:
    ok("era raw mean %d-%d" % (a, b), S.loc[a:b, "total"].mean(), v, .006)

print("\n--- Section 3.2 and Table 1, trends with no-trend bands ---")
for c, a, b, sl, R, lo_, hi_ in [
        ("total", 1884, 2023, 0.0584, 0.2036, -0.050, 0.050),
        ("any_int", 1884, 2023, 0.2042, 0.7197, -0.113, 0.115),
        ("pos_only", 1884, 2023, -0.1458, 0.6186, -0.087, 0.088),
        ("total", 1951, 2023, -0.1086, 0.2366, -0.090, 0.086),
        ("any_int", 1951, 2023, 0.0131, 0.0060, -0.037, 0.040),
        ("pos_only", 1951, 2023, -0.1217, 0.5513, -0.092, 0.087),
        ("total", 1977, 2023, -0.1026, 0.1474, -0.099, 0.099),
        ("any_int", 1977, 2023, -0.0192, 0.0058, -0.073, 0.078),
        ("pos_only", 1977, 2023, -0.0834, 0.4439, -0.081, 0.077)]:
    s = S.loc[a:b, c]
    ok("slope %s %d-%d" % (c, a, b), slope(s), sl, .0005)
    ok("R2 %s %d-%d" % (c, a, b), r2(s), R, .001)
    blo, bhi = band(s)
    ok("band low %s %d-%d" % (c, a, b), blo, lo_, .004)
    ok("band high %s %d-%d" % (c, a, b), bhi, hi_, .004)

print("\n--- Section 3.2, start-year sweep ---")
ok("start years sampled", len(sw), 59, 0)
ok("significant increases", (sw.sig == "+").sum(), 5, 0)
ok("significant decreases", (sw.sig == "-").sum(), 15, 0)
ok("minimum slope", sw.slope.min(), -0.1455, .0005)
ok("maximum slope", sw.slope.max(), 0.0584, .0005)
ok("last increasing start year", sw[sw.sig == "+"].start.max(), 1892, 0)
ok("first decreasing start year", sw[sw.sig == "-"].start.min(), 1946, 0)
ok("last decreasing start year", sw[sw.sig == "-"].start.max(), 1978, 0)
ok("classified start years", sw.csig.notna().sum(), 25, 0)
ok("classified significant", sw.csig.dropna().astype(str).isin(["+", "-"]).sum(), 0, 0)
sub = sw[sw.start <= 1994]
up = int((sub.sig == "+").sum()); dn = int((sub.sig == "-").sum())
ok("30-yr-minimum sweep, increases", up, 5, 0)
ok("30-yr-minimum sweep, decreases", dn, 15, 0)

print("\n--- Section 3.2, detection threshold and power ---")
s = S.loc[1951:2023, "any_int"]
x = np.asarray(s.index, float); y = np.asarray(s, float)
lo, hi = band(s)
resid = y - np.poly1d(np.polyfit(x, y, 1))(x); xc = x - x.mean()
for beta, want, tol in [(0.06, 84, 1), (0.08, 97, 1)]:
    hits = sum(1 for _ in range(600)
               if not (lo <= st.linregress(x, resid[rng.permutation(len(resid))]
                                           + beta * xc + y.mean()).slope <= hi))
    ok("power at %.2f storms/yr (%%)" % beta, 100 * hits / 600, want, tol)
ok("threshold expressed over 73 years", 0.06 * 73, 4.4, .05)
ok("observed slope as fraction of threshold", slope(s) / 0.06, 0.22, .02)
ok("raw decline as multiple of threshold", abs(slope(S.loc[1951:2023, "total"])) / 0.06, 1.8, .05)
print("\n--- Section 2, block-length sensitivity ---")


def blocklen_checks():
    windows = [("total", 1884, 2023), ("total", 1951, 2023),
               ("any_int", 1951, 2023), ("pos_only", 1951, 2023)]
    for L in (5, 15, 20):
        for cls, a, b in windows:
            ser = S.loc[a:b, cls]
            sl = slope(ser)
            lo10, hi10 = band(ser, L=10)
            loL, hiL = band(ser, L=L)
            same = (not (lo10 <= sl <= hi10)) == (not (loL <= sl <= hiL))
            # documented exception, see Section 2 of the manuscript
            want = 0 if (L == 20 and cls == "total" and a == 1884) else 1
            ok("L=%d verdict %s %d-%d matches L=10" % (L, cls, a, b),
               1.0 if same else 0.0, want, 0)


blocklen_checks()

print("\n--- Section 3.3, the identity closes ---")
ok("classified + unclassified, 1884-2023",
   slope(S.any_int) + slope(S.pos_only), 0.0584, .0005)
ok("classified + unclassified, 1951-2023",
   slope(S.loc[1951:2023, "any_int"]) + slope(S.loc[1951:2023, "pos_only"]), -0.1086, .0005)
for a, b, v in [(1945, 1950, 14.7), (1951, 1976, 18.6), (1977, 2000, 20.0), (2001, 2023, 18.9)]:
    ok("era classified %d-%d" % (a, b), S.loc[a:b, "any_int"].mean(), v, .05)
for a, b, v in [(1945, 1950, 6.3), (1951, 1976, 6.5), (1977, 2000, 2.1), (2001, 2023, 0.4)]:
    ok("era unclassified %d-%d" % (a, b), S.loc[a:b, "pos_only"].mean(), v, .05)

print("\n--- Section 3.4 and Table 2, backtest at both horizons ---")
for H, want in [(20, {"Full-record mean": (0.95, 1.15, -0.80),
                      "Trailing mean": (1.09, 1.18, -0.36),
                      "Linear trend extrapolated": (1.38, 1.69, 0.64),
                      "Last observed value": (2.79, 3.55, 0.02)}),
                (30, {"Full-record mean": (0.95, 0.99, -0.95),
                      "Trailing mean": (0.75, 0.79, -0.69),
                      "Linear trend extrapolated": (0.53, 0.82, 0.30),
                      "Last observed value": (2.13, 3.03, 0.55)})]:
    org = [t for t in range(1951 + H, 2024 - H)]
    ok("origins at %d-year horizon" % H, len(org), 33 if H == 20 else 13, 0)
    E = {k: [] for k in want}
    for t in org:
        tr = s[s.index <= t]; tg = s[(s.index > t) & (s.index <= t + H)].mean()
        lr = st.linregress(np.asarray(tr.index, float), np.asarray(tr, float))
        E["Full-record mean"].append(tr.mean() - tg)
        E["Trailing mean"].append(tr.iloc[-H:].mean() - tg)
        E["Linear trend extrapolated"].append(lr.slope * (t + (H + 1) / 2) + lr.intercept - tg)
        E["Last observed value"].append(tr.iloc[-1] - tg)
    for k, (mae, rmse, bias) in want.items():
        v = np.array(E[k])
        ok("H=%d %s MAE" % (H, k), np.abs(v).mean(), mae, .006)
        ok("H=%d %s RMSE" % (H, k), np.sqrt((v ** 2).mean()), rmse, .006)
        ok("H=%d %s bias" % (H, k), v.mean(), bias, .006)

print("\n--- Section 3.4, why neither horizon is decisive ---")
for H, want in [(20, 2.6), (30, 1.4)]:
    org = [t for t in range(1951 + H, 2024 - H)]
    ok("independent windows at H=%d" % H, ((org[-1] + H) - (org[0] + 1) + 1) / H, want, .06)
for lo_y, hi_y, wclim, wlin in [(1981, 1987, 1.00, 0.22), (1988, 1993, 0.88, 0.88)]:
    ec, el = [], []
    for t in range(lo_y, hi_y + 1):
        tr = s[s.index <= t]; tg = s[(s.index > t) & (s.index <= t + 30)].mean()
        lr = st.linregress(np.asarray(tr.index, float), np.asarray(tr, float))
        ec.append(tr.mean() - tg); el.append(lr.slope * (t + 15.5) + lr.intercept - tg)
    ok("origins %d-%d, climatology MAE" % (lo_y, hi_y), np.abs(np.array(ec)).mean(), wclim, .006)
    ok("origins %d-%d, linear trend MAE" % (lo_y, hi_y), np.abs(np.array(el)).mean(), wlin, .006)
ok("disjoint block mean 1951-1980", s.loc[1951:1980].mean(), 18.60, .006)
ok("disjoint block mean 1981-2010", s.loc[1981:2010].mean(), 19.50, .006)

print("\n--- Section 3.4, the level ---")
m = S.loc[1977:2023, "any_int"].mean()
ok("classified mean 1977-2023", m, 19.45, .006)
plo, phi = st.poisson.ppf([.025, .975], m)
ok("Poisson interval low", plo, 11, 0); ok("Poisson interval high", phi, 29, 0)

if CLIP:
    print("\n--- Table 3, sensitivity to the counting rule ---")
    d = pd.read_csv(CLIP, low_memory=False)
    d = d[(d.SEASON >= 1884) & (d.SEASON <= 2023)].copy()
    gcol = "TOKYO_GRADE" if "TOKYO_GRADE" in d.columns else "TOK_GRADE"
    grade = pd.to_numeric(d[gcol], errors="coerce").fillna(0)
    w = pd.to_numeric(d["USA_WIND"], errors="coerce").fillna(0)
    d["_any"] = grade.isin([2, 3, 4, 5, 9]) | (w > 0)
    main = d[d.TRACK_TYPE == "main"]
    cnt = main.groupby("SID").size()
    trop = main.groupby("SID").NATURE.agg(lambda v: bool(set(v) & {"TS", "DS"}))
    Y = range(1884, 2024)
    rules = {
        "Main track only (primary)": (main, 0.0584, -0.1086, 0.0131),
        "All SIDs including spurs": (d, 0.0608, -0.1215, 0.0108),
        "Main track, tropical nature": (main[main.SID.isin(trop[trop].index)], 0.0856, -0.1128, 0.0089),
        "Main track, 2+ PAR fixes": (main[main.SID.isin(cnt[cnt >= 2].index)], 0.0569, -0.1060, 0.0092)}
    for lab, (sub, a1, a2, a3) in rules.items():
        n = lambda x: x.groupby("SEASON").SID.nunique().reindex(Y).fillna(0).astype(int)
        tot, cls = n(sub), n(sub[sub._any])
        ok("%s: raw 1884-2023" % lab, slope(tot), a1, .0005)
        ok("%s: raw 1951-2023" % lab, slope(tot.loc[1951:2023]), a2, .0005)
        ok("%s: classified 1951-2023" % lab, slope(cls.loc[1951:2023]), a3, .0005)
else:
    print("\n--- Table 3 skipped: pass par_clipped.csv as the third argument to enable ---")

print("\n%d checks run, %d failed." % (N[0], FAIL[0]))
sys.exit(1 if FAIL[0] else 0)
