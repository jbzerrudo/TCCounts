"""
Version 4, step 3. When the record changed.

Two families of series, main track, 1951 to 2023 (HKO from 1961):

  shares  each agency's annual share of the storms in a reference set, for all
          systems (reference: every main-track storm) and at tropical-storm
          strength (reference: storms any agency rates tropical storm or
          stronger); binomial
  counts  annual counts of each storm class from decompose.py; Poisson

Each series is fitted with a constant, a smooth trend (logit-linear for shares,
log-linear for counts), and piecewise-constant means with 1, 2 or 3 breaks
(every segment at least 5 years, exhaustive search). Models are compared by a
quasi-BIC that divides the log-likelihood by the overdispersion of the 3-break
fit (floored at 1). For the selected step model, each break's support interval
holds the other breaks fixed and keeps positions whose scaled deviance is
within 3.84 of the best.

Outputs
    data/dating_shares.csv, data/dating_counts.csv

Usage
    python src/dating.py data/storm_flags.csv data/
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

from classes import annual, any4, load_flags
from decompose import classes

MINSEG = 5
COUNT_CLASSES = ("ts_or_stronger", "depression_only", "no_intensity", "cma_hko_only",
                 "jtwc_not_jma", "ts_all3")


def _xlogy(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(a > 0, a * np.log(np.where(b > 0, b, 1.0)), 0.0)


def _segments(L, m):
    for br in itertools.combinations(range(MINSEG, L - MINSEG + 1), m):
        if all(b2 - b1 >= MINSEG for b1, b2 in zip(br[:-1], br[1:])):
            yield br


def _fill(L, br, values):
    edges = [0, *br, L]
    out = np.empty(L)
    for (a, b), v in zip(zip(edges[:-1], edges[1:]), values):
        out[a:b] = v
    return out


class Binomial:
    def __init__(self, k, n):
        self.k, self.n = k, n

    def ll(self, p):
        return float((_xlogy(self.k, p) + _xlogy(self.n - self.k, 1 - p)).sum())

    def step(self, br):
        edges = [0, *br, len(self.k)]
        vals = [self.k[a:b].sum() / self.n[a:b].sum() for a, b in zip(edges[:-1], edges[1:])]
        return _fill(len(self.k), br, vals)

    def smooth(self, t):
        x = (t - t.mean()) / 10.0
        X = np.c_[np.ones_like(x), x]
        p0 = np.clip(self.k.sum() / self.n.sum(), 1e-6, 1 - 1e-6)
        beta = np.array([np.log(p0 / (1 - p0)), 0.0])
        for _ in range(100):
            p = 1 / (1 + np.exp(-(X @ beta)))
            W = self.n * p * (1 - p)
            z = X @ beta + (self.k - self.n * p) / np.maximum(W, 1e-12)
            new = np.linalg.solve(X.T @ (W[:, None] * X), X.T @ (W * z))
            done = np.max(np.abs(new - beta)) < 1e-10
            beta = new
            if done:
                break
        return 1 / (1 + np.exp(-(X @ beta))), beta[1]

    def pearson(self, p):
        v = self.n * p * (1 - p)
        ok = v > 0
        return float((((self.k - self.n * p) ** 2)[ok] / v[ok]).sum())

    def label(self, a, b):
        return f"{100 * self.k[a:b].sum() / self.n[a:b].sum():.1f}%"


class Poisson:
    def __init__(self, k):
        self.k = k

    def ll(self, mu):
        return float(_xlogy(self.k, mu).sum() - mu.sum())

    def step(self, br):
        edges = [0, *br, len(self.k)]
        return _fill(len(self.k), br, [self.k[a:b].mean() for a, b in zip(edges[:-1], edges[1:])])

    def smooth(self, t):
        x = (t - t.mean()) / 10.0
        X = np.c_[np.ones_like(x), x]
        beta = np.array([np.log(max(self.k.mean(), 1e-6)), 0.0])
        for _ in range(100):
            mu = np.exp(X @ beta)
            z = X @ beta + (self.k - mu) / mu
            new = np.linalg.solve(X.T @ (mu[:, None] * X), X.T @ (mu * z))
            done = np.max(np.abs(new - beta)) < 1e-10
            beta = new
            if done:
                break
        return np.exp(X @ beta), beta[1]

    def pearson(self, mu):
        ok = mu > 0
        return float((((self.k - mu) ** 2)[ok] / mu[ok]).sum())

    def label(self, a, b):
        return f"{self.k[a:b].mean():.2f}/yr"


def date_series(model, t):
    L = len(t)
    fits = {"constant": (model.ll(model.step(())), 1, ())}
    p_s, slope10 = model.smooth(t)
    fits["smooth"] = (model.ll(p_s), 2, ())
    best_fit = {}
    for m in (1, 2, 3):
        top = max(((model.ll(model.step(br)), br) for br in _segments(L, m)), key=lambda v: v[0])
        fits[f"{m} break"] = (top[0], 2 * m + 1, top[1])
        best_fit[m] = model.step(top[1])
    phi = max(1.0, model.pearson(best_fit[3]) / (L - 7))
    q = {k: -2 * f[0] / phi + f[1] * np.log(L) for k, f in fits.items()}
    sel = min(q, key=q.get)
    llb, _, br = fits[sel]
    support = []
    for i in range(len(br)):
        others = list(br[:i]) + list(br[i + 1:])
        ok = []
        for pos in range(MINSEG, L - MINSEG + 1):
            cand = sorted(others + [pos])
            if any(b2 - b1 < MINSEG for b1, b2 in zip([0] + cand, cand + [L])):
                continue
            if 2 * (llb - model.ll(model.step(tuple(cand)))) / phi <= 3.84:
                ok.append(pos)
        support.append(f"{t[min(ok)]}-{t[max(ok)]}")
    edges = [0, *br, L]
    segs = "; ".join(f"{t[a]}-{t[b - 1]}: {model.label(a, b)}" for a, b in zip(edges[:-1], edges[1:]))
    best_step = min(q[f"{m} break"] for m in (1, 2, 3))
    return dict(selected=sel, breaks=" ".join(str(t[b]) for b in br), support=" ".join(support),
                segments=segs, phi=round(phi, 2),
                dqbic_smooth_minus_selected=round(q["smooth"] - q[sel], 1),
                dqbic_smooth_minus_best_step=round(q["smooth"] - best_step, 1),
                smooth_slope_per_decade=round(slope10, 3))


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "data/storm_flags.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "data"
    shares, counts = [], []
    for dom in ("WNP basin", "PAR"):
        fl = load_flags(src, dom)
        ts = any4(fl, "TS")
        for cls in ("all systems", "TS strength"):
            ref = pd.Series(True, index=fl.index) if cls == "all systems" else ts
            for ag in ("JMA", "JTWC", "CMA", "HKO"):
                col = f"{ag}_all" if cls == "all systems" else f"{ag}_TS"
                y0 = 1961 if ag == "HKO" else 1951
                t = np.arange(y0, 2024)
                n = annual(fl, ref).loc[y0:2023].to_numpy(float)
                k = annual(fl, ref & fl[col]).loc[y0:2023].to_numpy(float)
                r = date_series(Binomial(k, n), t)
                shares.append(dict(domain=dom, cls=cls, agency=ag, years=f"{y0}-2023", **r))
        cl = classes(fl)
        t = np.arange(1951, 2024)
        for name in COUNT_CLASSES:
            k = annual(fl, cl[name]).loc[1951:2023].to_numpy(float)
            counts.append(dict(domain=dom, cls=name, **date_series(Poisson(k), t)))
    S, C = pd.DataFrame(shares), pd.DataFrame(counts)
    S.to_csv(os.path.join(out, "dating_shares.csv"), index=False)
    C.to_csv(os.path.join(out, "dating_counts.csv"), index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 90)
    print(C[["domain", "cls", "selected", "segments"]].to_string(index=False))
    print(f"written: {out}/dating_shares.csv, {out}/dating_counts.csv")


if __name__ == "__main__":
    main()
