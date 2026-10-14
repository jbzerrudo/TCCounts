"""
Check every value reported from the version 4 analyses against the archived
outputs. Prints OK or MISMATCH beside each value and exits non-zero if any
check fails, like analysis.py.

Deterministic values are checked to the rounding they were reported at: a
value reported as -0.053 must lie within 0.0005 of it. Simulated values (test
size, power, calibrated p-values, the curve-level test) carry Monte Carlo
tolerances, so the checks also pass on a numpy or statsmodels version whose
random streams or optimiser differ from the one used here.

Usage
    python src/check_v4.py data/
"""
import os
import sys

import numpy as np
import pandas as pd

from classes import REQUIREMENTS, annual, fit_many, load_flags, requirement

D = sys.argv[1] if len(sys.argv) > 1 else "data"
rd = lambda f: pd.read_csv(os.path.join(D, f))
N, FAIL = [0], [0]


def _report(name, good, want, got):
    N[0] += 1
    FAIL[0] += 0 if good else 1
    print(("  OK      " if good else "  MISMATCH") + f"  {name:<70s} reported={want:<12} computed={got}")


def ok(name, got, reported):
    """reported is the value as written (a string), so its rounding sets the tolerance."""
    dec = len(reported.split(".")[1]) if "." in reported else 0
    tol = 0.5 * 10 ** -dec + 1e-9
    _report(name, abs(float(got) - float(reported)) <= tol, reported, round(float(got), dec + 2))


def exact(name, got, want):
    _report(name, got == want, str(want), got)


def within(name, got, lo, hi):
    _report(name, lo <= float(got) <= hi, f"{lo:g}..{hi:g}", round(float(got), 4))


G = rd("grid_threshold.csv")
V = {"-": -1, "0": 0, "+": 1}
pct = lambda g: 100 * (g.verdict == "-").mean()
cell = lambda dom, req, t: G[(G.domain == dom) & (G.requirement == req) & (G.threshold == t)]
main = lambda dom, req, t, y: G[(G.domain == dom) & (G.requirement == req) & (G.threshold == t) &
                                (G.rule == "main track") & (G.start == y)].iloc[0]

print("\n--- Consistency with the published grid and the published test ---")
M = rd("multiverse.csv")
for new, old in [("JMA", "JMA"), ("JTWC", "JTWC"), ("JMA|JTWC", "JMA+JTWC"),
                 ("any of 3", "JMA+JTWC+CMA"), ("any of 4", "all four"), ("raw total", "raw total")]:
    a = G[(G.requirement == new) & (G.threshold.isin(["all", "-"]))].sort_values(["domain", "rule", "start"])
    b = M[M.agency == old].sort_values(["domain", "rule", "start"])
    within(f"{new}: largest slope difference from multiverse.csv",
           np.abs(a.slope.values - b.slope.values).max(), 0, 1e-12)
    exact(f"{new}: verdicts that differ from multiverse.csv", int((a.verdict.values != b.verdict.values).sum()), 0)

mism, n_cmp = 0, 0
for dom in ("WNP basin", "PAR"):
    for rule in ("all SIDs", "main track", ">=2 fixes", "tropical nature"):
        fl = load_flags(os.path.join(D, "storm_flags.csv"), dom, rule)
        names = [(r, t) for t in ("all", "TD", "TS") for r in REQUIREMENTS] + [("raw total", "-")]
        Mx = np.vstack([annual(fl, pd.Series(True, index=fl.index)).values if r == "raw total"
                        else annual(fl, requirement(fl, r, t)).values for r, t in names])
        sub = G[(G.domain == dom) & (G.rule == rule)].set_index(["requirement", "threshold", "start"])
        for y in range(1951, 2001):
            _, _, vd = fit_many(Mx[:, y - 1884:], np.arange(y, 2024))
            for (r, t), v in zip(names, vd):
                n_cmp += 1
                mism += int(V[sub.loc[(r, t, y), "verdict"]] != v)
exact("specifications compared between fit_many and multiverse.fit", n_cmp, 7600)
exact("verdicts that differ between fit_many and multiverse.fit", mism, 0)

print("\n--- Archive ---")
A = rd("archive_meta.csv").set_index("TRACK_TYPE")
exact("PROVISIONAL rows in 1884-2023", int(A.loc["PROVISIONAL", "n_1884_2023"]), 0)
exact("US-PROVISIONAL rows in 1884-2023", int(A.loc["US-PROVISIONAL", "n_1884_2023"]), 0)
exact("PROVISIONAL and US-PROVISIONAL rows exist after 2023",
      bool(A.loc["PROVISIONAL", "n_after_2023"] > 0 and A.loc["US-PROVISIONAL", "n_after_2023"] > 0), True)
C = rd("cma_consistency.csv").iloc[0]
ok("CMA fixes whose category disagrees with their wind (%)",
   100 * (C.ts_category_wind_below_33kt + C.td_or_weaker_category_wind_34kt_plus) / C.fixes_with_cma_wind, "1")

print("\n--- Analysis 1: share of specifications with a significant decline, and median slope ---")
T1 = {("WNP basin", "all"): [("27.5", "-0.086"), ("11.5", "-0.005"), ("75.5", "-0.173"), ("12.5", "-0.019"), ("87.0", "-0.174"), ("16.0", "-0.073")],
      ("WNP basin", "TS"): [("27.0", "-0.085"), ("12.0", "-0.070"), ("41.5", "-0.100"), ("13.5", "-0.066"), ("44.0", "-0.088"), ("20.0", "-0.085")],
      ("PAR", "all"): [("12.0", "-0.041"), ("0.0", "-0.007"), ("54.5", "-0.104"), ("0.0", "-0.014"), ("62.5", "-0.101"), ("0.0", "-0.025")],
      ("PAR", "TS"): [("0.0", "-0.034"), ("19.0", "-0.051"), ("28.0", "-0.055"), ("0.0", "-0.030"), ("19.5", "-0.044"), ("28.5", "-0.056")]}
for (dom, t), vals in T1.items():
    for req, (p, m) in zip(REQUIREMENTS, vals):
        g = cell(dom, req, t)
        ok(f"{dom}, {req}, {t}: % significant declines", pct(g), p)
        ok(f"{dom}, {req}, {t}: median slope", g.slope.median(), m)
exact("PAR, JTWC, all systems: significant increases", int((cell("PAR", "JTWC", "all").verdict == "+").sum()), 1)
for dom, s_all, s_ts in [("WNP basin", "0.169", "0.034"), ("PAR", "0.097", "0.027")]:
    for t, want in (("all", s_all), ("TS", s_ts)):
        med = [cell(dom, r, t).slope.median() for r in REQUIREMENTS]
        ok(f"{dom}, {t}: spread of median slopes across requirements", max(med) - min(med), want)
for dom, n_all, n_ts in [("WNP basin", "81.9", "93.8"), ("PAR", "78.8", "93.1")]:
    for t, want in (("all", n_all), ("TS", n_ts)):
        g = G[(G.domain == dom) & (G.threshold == t) & G.requirement.isin(REQUIREMENTS)]
        ok(f"{dom}, {t}: % of specifications with a negative slope", 100 * (g.slope < 0).mean(), want)


def eta2(g, f):
    return 100 * ((g.groupby(f).slope.transform("mean") - g.slope.mean()) ** 2).sum() / ((g.slope - g.slope.mean()) ** 2).sum()


for dom, ag_all, ag_ts, st_ts in [("WNP basin", "29", "2", "70"), ("PAR", "27", "15", "72")]:
    for t, want in (("all", ag_all), ("TS", ag_ts)):
        g = G[(G.domain == dom) & (G.threshold == t) & G.requirement.isin(REQUIREMENTS)]
        ok(f"{dom}, {t}: % of slope variance explained by the requirement", eta2(g, "requirement"), want)
    g = G[(G.domain == dom) & (G.threshold == "TS") & G.requirement.isin(REQUIREMENTS)]
    ok(f"{dom}, TS: % of slope variance explained by the start year", eta2(g, "start"), st_ts)

print("\n--- Analysis 2: where the 1951-2023 decline sits (main track) ---")
P = rd("decomposition.csv")
E = rd("era_means.csv")
dget = lambda dom, y, c: P[(P.domain == dom) & (P.start == y) & (P.cls == c)].iloc[0]
eget = lambda dom, c, era: E[(E.domain == dom) & (E.cls == c) & (E.era == era)].mean_per_year.iloc[0]
T2 = {"WNP basin": [("ts_or_stronger", "-0.053", "0", "32", "29.5", "26.1"),
                    ("depression_only", "-0.080", "-", "49", "8.4", "3.8"),
                    ("below_depression", "0.001", "0", "-1", "0.04", "0.09"),
                    ("no_intensity", "-0.033", "-", "20", "1.6", "0.0")],
      "PAR": [("ts_or_stronger", "-0.039", "-", "36", "17.1", "15.2"),
              ("depression_only", "-0.044", "-", "40", "6.2", "3.5"),
              ("below_depression", "0.006", "0", "-6", "0.35", "0.57"),
              ("no_intensity", "-0.0315", "-", "29", "1.5", "0.0")]}
for dom, rows in T2.items():
    parts = 0.0
    for c, sl, vd, sh, m0, m1 in rows:
        r = dget(dom, 1951, c)
        parts += r.slope
        ok(f"{dom}, {c}: slope", r.slope, sl)
        exact(f"{dom}, {c}: verdict", r.verdict, vd)
        ok(f"{dom}, {c}: share of the total slope (%)", r.share_of_total_pct, sh)
        ok(f"{dom}, {c}: storms per year 1951-1976", eget(dom, c, "1951-1976"), m0)
        ok(f"{dom}, {c}: storms per year 2001-2023", eget(dom, c, "2001-2023"), m1)
    within(f"{dom}: the four classes add up to the total slope (difference)",
           abs(parts - dget(dom, 1951, "total").slope), 0, 1e-12)
for dom, below, ts_drop, dep_drop in [("WNP basin", "68", "11", "55"), ("PAR", "64", "11", "43")]:
    g = P[(P.domain == dom) & (P.start == 1951)].set_index("cls")
    ok(f"{dom}: % of the decline below tropical-storm strength or ungraded",
       g.loc[["depression_only", "below_depression", "no_intensity"], "share_of_total_pct"].sum(), below)
    for c, want in (("ts_or_stronger", ts_drop), ("depression_only", dep_drop)):
        a, b = eget(dom, c, "1951-1976"), eget(dom, c, "2001-2023")
        ok(f"{dom}, {c}: % fall from 1951-1976 to 2001-2023", 100 * (a - b) / a, want)
for dom, want in (("WNP basin", "13.5"), ("PAR", "17.6")):
    r = dget(dom, 1951, "ts_or_stronger")
    ok(f"{dom}: tropical-storm decline from 1951 as % of the mean", -100 * r.slope * 72 / r.mean_per_year, want)

print("\n--- Analysis 1: tropical-storm-strength slopes by start year (main track) ---")
T3 = {"WNP basin": {1951: [("-0.043", "0"), ("-0.004", "0"), ("-0.073", "-"), ("-0.053", "0")],
                    1961: [("-0.098", "-"), ("-0.072", "0"), ("-0.132", "-"), ("-0.104", "-")],
                    1977: [("-0.080", "0"), ("-0.061", "0"), ("-0.081", "0"), ("-0.065", "0")]},
      "PAR": {1951: [("-0.023", "0"), ("-0.035", "0"), ("-0.048", "-"), ("-0.039", "-")],
              1961: [("-0.036", "0"), ("-0.049", "-"), ("-0.058", "-"), ("-0.045", "0")],
              1977: [("-0.036", "0"), ("-0.051", "0"), ("-0.053", "0"), ("-0.039", "0")]}}
for dom, by in T3.items():
    for y, vals in by.items():
        for req, (sl, vd) in zip(("JMA", "JTWC", "CMA", "any of 3"), vals):
            r = main(dom, req, "TS", y)
            ok(f"{dom}, {req}, TS, from {y}: slope", r.slope, sl)
            exact(f"{dom}, {req}, TS, from {y}: verdict", r.verdict, vd)
r = main("WNP basin", "HKO", "TS", 1961)
ok("WNP basin, HKO, TS, from 1961: slope", r.slope, "-0.115")
exact("WNP basin, HKO, TS, from 1961: verdict", r.verdict, "-")

print("\n--- Analysis 1: coverage at tropical-storm strength, basin, by era ---")
CV = rd("ts_coverage.csv")
for ag, vals in (("JMA", ("91.4", "93.1", "92.5")), ("JTWC", ("87.9", "95.5", "94.7")), ("CMA", ("97.1", "93.0", "93.2"))):
    for era, want in zip(("1951-1976", "1977-2000", "2001-2023"), vals):
        g = CV[(CV.domain == "WNP basin") & (CV.agency == ag) & (CV.era == era)]
        ok(f"WNP basin, {ag}: % of tropical-storm-strength storms rated TS, {era}", g.pct.iloc[0], want)

print("\n--- Analysis 1: other items ---")
r = main("PAR", "JTWC", "TD", 1951)
ok("PAR, JTWC, depression or stronger, from 1951: slope", r.slope, "0.018")
exact("PAR, JTWC, depression or stronger, from 1951: verdict", r.verdict, "0")
K = rd("class_counts.csv")
h = K[(K.domain == "WNP basin") & (K.SEASON >= 1951)]
exact("WNP basin: storms graded only by HKO, 1951-2023", int(h.hko_only.sum()), 24)
exact("WNP basin: storms graded only by HKO, 1960s", int(h[(h.SEASON >= 1960) & (h.SEASON <= 1969)].hko_only.sum()), 14)
for dom in ("WNP basin", "PAR"):
    for t in ("all", "TD", "TS"):
        exact(f"{dom}, all 3, {t}, main track, from 1951: verdict", main(dom, "all 3", t, 1951).verdict, "0")
W = rd("cma_wind_sensitivity.csv")
diffs = []
for dom in ("WNP basin", "PAR"):
    for a_, b_ in (("CMA, TS by category", "CMA, TS by wind"), ("any of 3, TS, CMA by category", "any of 3, TS, CMA by wind")):
        for y in (1951, 1961, 1977):
            ra = W[(W.domain == dom) & (W.series == a_) & (W.start == y)].iloc[0]
            rb = W[(W.domain == dom) & (W.series == b_) & (W.start == y)].iloc[0]
            diffs.append((abs(ra.slope - rb.slope), ra.verdict != rb.verdict, dom, a_, y))
ok("CMA by wind instead of category: largest change in slope", max(d[0] for d in diffs), "0.006")
flips = [d for d in diffs if d[1]]
exact("CMA by wind instead of category: verdicts that flip", len(flips), 1)
exact("CMA by wind instead of category: the flip is PAR, any of 3, from 1961",
      bool(len(flips) == 1 and flips[0][2] == "PAR" and flips[0][3].startswith("any of 3") and flips[0][4] == 1961), True)

print("\n--- Analysis 4: when the record changed ---")
DC = rd("dating_counts.csv").set_index(["domain", "cls"])
T4 = {("WNP basin", "ts_or_stronger"): ("", ["28.1"]),
      ("WNP basin", "depression_only"): ("1967 1972 2004", ["7.6", "14.0", "6.0", "3.6"]),
      ("WNP basin", "no_intensity"): ("1971", ["2.05", "0.04"]),
      ("WNP basin", "cma_hko_only"): ("1967 1972 1991", ["9.1", "16.0", "5.2", "0.8"]),
      ("WNP basin", "jtwc_not_jma"): ("1993 2004", ["1.6", "9.4", "4.8"]),
      ("PAR", "ts_or_stronger"): ("", ["16.2"]),
      ("PAR", "depression_only"): ("2003", ["5.9", "3.3"]),
      ("PAR", "no_intensity"): ("1972 1986", ["1.81", "0.29", "0.00"]),
      ("PAR", "cma_hko_only"): ("1990", ["4.5", "0.35"]),
      ("PAR", "jtwc_not_jma"): ("1990 2003", ["0.85", "4.85", "2.4"])}
for (dom, c), (brk, means) in T4.items():
    r = DC.loc[(dom, c)]
    got = "" if pd.isna(r.breaks) else str(r.breaks)
    exact(f"{dom}, {c}: breaks", got if got else r.selected, brk if brk else "constant")
    segs = [float(p.split(": ")[1].replace("/yr", "")) for p in r.segments.split("; ")]
    exact(f"{dom}, {c}: number of segments", len(segs), len(means))
    for i, (g_, w_) in enumerate(zip(segs, means)):
        ok(f"{dom}, {c}: segment {i + 1}, storms per year", g_, w_)
DS = rd("dating_shares.csv").set_index(["domain", "cls", "agency"])
share_segs = lambda r: [float(p.split(": ")[1].rstrip("%")) for p in r.segments.split("; ")]
r = DS.loc[("WNP basin", "all systems", "JTWC")]
exact("WNP basin, JTWC share of all systems: breaks", str(r.breaks), "1972 1991")
for i, (g_, w_) in enumerate(zip(share_segs(r), ("63.7", "82.6", "96.7"))):
    ok(f"WNP basin, JTWC share of all systems: segment {i + 1} (%)", g_, w_)
ok("WNP basin, JTWC share: smooth trend worse by (QBIC)", r.dqbic_smooth_minus_selected, "38")
r = DS.loc[("WNP basin", "all systems", "CMA")]
exact("WNP basin, CMA share of all systems: breaks", str(r.breaks), "1967 1993 2004")
for i, (g_, w_) in enumerate(zip(share_segs(r), ("90.4", "97.2", "81.3", "93.2"))):
    ok(f"WNP basin, CMA share of all systems: segment {i + 1} (%)", g_, w_)
ok("WNP basin, CMA share: smooth trend worse by (QBIC)", r.dqbic_smooth_minus_selected, "45")
r = DS.loc[("PAR", "TS strength", "JTWC")]
exact("PAR, JTWC share at tropical-storm strength: model", r.selected, "constant")
ok("PAR, JTWC share at tropical-storm strength (%)", share_segs(r)[0], "91.9")

print("\n--- Analysis 5a: size of the test (simulated; Monte Carlo tolerance) ---")
S1 = rd("size_ar1.csv")
s = lambda dom, ser, n: S1[(S1.domain == dom) & (S1.series == ser) & (S1.n == n)].size_pct.iloc[0]
within("WNP basin total, AR(1), 73 yr: size (%), reported 6.8", s("WNP basin", "total", 73), 5.3, 8.3)
within("WNP basin total, AR(1), 24 yr: size (%), reported 9.6", s("WNP basin", "total", 24), 8.1, 11.1)
p_ = S1[S1.domain == "PAR"].size_pct
within("PAR, AR(1), both series, all windows: smallest size (%), reported 3.7", p_.min(), 2.7, 4.7)
within("PAR, AR(1), both series, all windows: largest size (%), reported 5.4", p_.max(), 4.4, 6.4)
SF = rd("size_fitted.csv")
f = SF[(SF.domain == "WNP basin") & (SF.series == "total")]
exact("WNP basin total: AIC chooses AR(3) with 0.28, 0.02, 0.25",
      bool(((f.model == "AR(3) by AIC") & (f.params == "ar=0.28,0.02,0.25")).any()), True)
fa = f[f.model == "AR(3) by AIC"].size_pct
within("WNP basin total, AR(3): smallest size over windows (%), reported 19", fa.min(), 17.5, 20.5)
within("WNP basin total, AR(3): largest size over windows (%), reported 21", fa.max(), 19.5, 22.5)
fb = f[f.model == "ARMA(1,1)"].size_pct
within("WNP basin total, ARMA(1,1): smallest size (%), reported about the same", fb.min(), 17.0, 20.5)
within("WNP basin total, ARMA(1,1): largest size (%), reported about the same", fb.max(), 19.5, 22.5)
f = SF[(SF.domain == "WNP basin") & (SF.series == "TS or stronger")]
exact("WNP basin TS: AIC chooses white noise", bool((f.model == "AR(0) by AIC").any()), True)
within("WNP basin TS, white noise: mean size (%), reported 4", f[f.model == "AR(0) by AIC"].size_pct.mean(), 3.0, 5.0)
fm = f[f.model == "ARMA(1,1)"].size_pct
within("WNP basin TS, ARMA(1,1): smallest size (%), reported 11", fm.min(), 9.5, 12.5)
within("WNP basin TS, ARMA(1,1): largest size (%), reported 13", fm.max(), 11.5, 14.5)
for ser in ("total", "TS or stronger"):
    f = SF[(SF.domain == "PAR") & (SF.series == ser)]
    exact(f"PAR {ser}: AIC chooses white noise", bool((f.model == "AR(0) by AIC").any()), True)
    within(f"PAR {ser}, white noise: mean size (%), reported 4", f[f.model == "AR(0) by AIC"].size_pct.mean(), 3.0, 5.0)

print("\n--- Analysis 5a: calibrated p-values (simulated; Monte Carlo tolerance) ---")
CP = rd("calibrated_p.csv")
T5 = [("WNP basin", "raw total", 1951, 0.008, 0.008), ("WNP basin", "raw total", 1961, 0.001, 0.001),
      ("WNP basin", "raw total", 1977, 0.23, 0.23),
      ("WNP basin", "TS or stronger", 1951, 0.046, 0.16), ("WNP basin", "TS or stronger", 1961, 0.002, 0.026),
      ("WNP basin", "TS or stronger", 1977, 0.20, 0.36),
      ("PAR", "raw total", 1951, 0.0, 0.001), ("PAR", "raw total", 1961, 0.0, 0.001),
      ("PAR", "raw total", 1977, 0.02, 0.03),
      ("PAR", "TS or stronger", 1951, 0.03, 0.03), ("PAR", "TS or stronger", 1961, 0.04, 0.05),
      ("PAR", "TS or stronger", 1977, 0.27, 0.27)]
for dom, ser, y, lo, hi in T5:
    g = CP[(CP.domain == dom) & (CP.series == ser) & (CP.start == y)].p_two_sided
    tol = 0.0015 if hi < 0.01 else max(0.006, 0.15 * hi)
    within(f"{dom}, {ser}, from {y}: smaller calibrated p (reported {lo:g}-{hi:g})", g.min(), max(0.0, lo - tol), hi + tol)
    within(f"{dom}, {ser}, from {y}: larger calibrated p (reported {lo:g}-{hi:g})", g.max(), max(0.0, lo - tol), hi + tol)
r = CP[(CP.domain == "WNP basin") & (CP.series == "raw total") & (CP.start == 1977)].iloc[0]
exact("WNP basin total from 1977: the lag-one test says significant decline", r.test_verdict, "-")

print("\n--- Analysis 5a: trend reaching 80% power (simulated; 4% relative tolerance) ---")
PW = rd("power.csv")
T6 = {("WNP basin", "total"): [(73, 0.141, 29), (53, 0.229, 36), (33, 0.470, 47)],
      ("WNP basin", "TS or stronger"): [(73, 0.091, 23), (53, 0.149, 28), (33, 0.304, 36)],
      ("PAR", "total"): [(73, 0.071, 23), (53, 0.116, 29), (33, 0.238, 38)],
      ("PAR", "TS or stronger"): [(73, 0.053, 24), (53, 0.087, 29), (33, 0.180, 37)]}
for (dom, ser), rows in T6.items():
    for n, b, chg in rows:
        r = PW[(PW.domain == dom) & (PW.series == ser) & (PW.n == n)].iloc[0]
        within(f"{dom}, {ser}, {n} yr: trend at 80% power (reported {b})", r.b80, 0.96 * b - 0.0005, 1.04 * b + 0.0005)
        within(f"{dom}, {ser}, {n} yr: change as % of window mean (reported {chg})",
               r.change_pct_of_mean, 0.96 * chg - 0.5, 1.04 * chg + 0.5)

print("\n--- Analysis 5b: curve-level test (2,000 synthetic archives; Monte Carlo tolerance) ---")
CT = rd("curve_test.csv")
T7 = {("WNP basin", "original"): (620, 314, 458, 0.005, 0.011), ("PAR", "original"): (400, 216, 371, 0.013, 0.043),
      ("WNP basin", "all"): (460, 300, 446, 0.014, 0.046), ("PAR", "all"): (258, 191, 268, 0.032, 0.055),
      ("WNP basin", "TS"): (316, 333, 444, 0.057, 0.096), ("PAR", "TS"): (190, 130, 195, 0.030, 0.052)}
for (dom, fam), (obs, q_lo, q_hi, p_lo, p_hi) in T7.items():
    g = CT[(CT.domain == dom) & (CT.family == fam)]
    exact(f"{dom}, {fam}: observed significant declines", int(g.obs_declines.iloc[0]), obs)
    within(f"{dom}, {fam}: null 95th percentile, smallest over block lengths (reported {q_lo})",
           g.null_declines_p95.min(), 0.85 * q_lo, 1.15 * q_lo)
    within(f"{dom}, {fam}: null 95th percentile, largest over block lengths (reported {q_hi})",
           g.null_declines_p95.max(), 0.85 * q_hi, 1.15 * q_hi)
    within(f"{dom}, {fam}: p, smallest over block lengths (reported {p_lo})", g.p_declines.min(),
           max(0, p_lo - 0.015), p_lo + 0.015)
    within(f"{dom}, {fam}: p, largest over block lengths (reported {p_hi})", g.p_declines.max(),
           max(0, p_hi - 0.02), p_hi + 0.02)
g = CT[(CT.domain == "WNP basin") & (CT.family == "original")]
within("WNP basin, original: null 95th percentile as % of 1,200, smallest (reported 26)",
       100 * g.null_declines_p95.min() / 1200, 22, 30)
within("WNP basin, original: null 95th percentile as % of 1,200, largest (reported 38)",
       100 * g.null_declines_p95.max() / 1200, 33, 44)

print("\n--- Analysis 5c: end year ---")
EY = rd("end_year.csv")
e = lambda dom, ser, end: EY[(EY.domain == dom) & (EY.series == ser) & (EY.end == end)].iloc[0]
for end, vd in ((2005, "0"), (2010, "0"), (2015, "-"), (2020, "-"), (2023, "-")):
    exact(f"WNP basin, any of 3, all systems, 1951-{end}: verdict",
          e("WNP basin", "any of 3, all systems", end).verdict_from_1951, vd)
r = e("WNP basin", "JMA|JTWC, all systems", 2005)
ok("WNP basin, JMA|JTWC, all systems, 1951-2005: slope", r.slope_from_1951, "0.135")
exact("WNP basin, JMA|JTWC, all systems, 1951-2005: verdict", r.verdict_from_1951, "+")
for end in (2005, 2010, 2015, 2020, 2023):
    exact(f"WNP basin, any of 3, TS, 1951-{end}: verdict", e("WNP basin", "any of 3, TS", end).verdict_from_1951, "0")
    exact(f"PAR, any of 3, TS, 1951-{end}: verdict", e("PAR", "any of 3, TS", end).verdict_from_1951,
          "-" if end == 2023 else "0")

print("\n%d checks run, %d failed." % (N[0], FAIL[0]))
sys.exit(1 if FAIL[0] else 0)
