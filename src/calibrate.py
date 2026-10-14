"""
Version 4, step 4. Size and power of the test, recalibrated.

Series: the main-track total and the count of storms any agency rates tropical
storm or stronger, in both domains; for the calibrated p-values also the
any-of-3 and JMA|JTWC counts (all-systems definition). Noise parameters come
from the 1951-2023 series after removing its OLS trend, so the trend no longer
inflates the dispersion, and the change a trend implies is expressed against
the mean of the window it is fitted to.

  size_ar1      rejection rate with no trend, Gaussian AR(1) noise with the
                detrended variance and lag-one autocorrelation
  power         trend reaching 80% power under that AR(1), windows of 24-73 yr
  size_fitted   rejection rate with no trend under the AR(p) chosen by AIC
                (p <= 5) and under an ARMA(1,1), both fitted by maximum
                likelihood to the detrended residuals
  calibrated_p  two-sided p-value of the observed slope under each fitted
                model, for windows starting 1951, 1961 and 1977

Requires statsmodels. Seeds are fixed; on another numpy version the simulated
values can differ within Monte Carlo error, which check_v4.py allows for.

Outputs
    data/size_ar1.csv, data/power.csv, data/size_fitted.csv, data/calibrated_p.csv

Usage
    python src/calibrate.py data/storm_flags.csv data/
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.signal import lfilter
from statsmodels.tsa.arima.model import ARIMA

from classes import annual, any4, fit_many, load_flags, requirement

warnings.filterwarnings("ignore")
NS = (24, 33, 43, 53, 63, 73)
R_SIZE, R_POW, R_P = 5000, 3000, 20000


def ar1(reps, n, sd, rho, rng):
    e = rng.standard_normal((reps, n))
    z = np.empty((reps, n))
    z[:, 0] = e[:, 0]
    s = np.sqrt(1 - rho ** 2)
    for i in range(1, n):
        z[:, i] = rho * z[:, i - 1] + s * e[:, i]
    return sd * z


def calib_rho(target, n=73, reps=4000, seed=7):
    """Latent AR(1) coefficient whose detrended lag-one autocorrelation matches target."""
    x = np.arange(n, dtype=float)

    def delivered(rho):
        Z = ar1(reps, n, 1.0, rho, np.random.default_rng(seed))
        sl, _, _ = fit_many(Z, x)
        R = Z - Z.mean(1, keepdims=True) - sl[:, None] * (x - x.mean())[None, :]
        a, c = R[:, :-1], R[:, 1:]
        am, cm = a.mean(1, keepdims=True), c.mean(1, keepdims=True)
        return float(np.mean(((a - am) * (c - cm)).sum(1) /
                             np.sqrt(((a - am) ** 2).sum(1) * ((c - cm) ** 2).sum(1))))
    if target <= 0:
        return 0.0
    lo, hi, best = 0.0, 0.98, None
    for _ in range(16):
        mid = 0.5 * (lo + hi)
        got = delivered(mid)
        if best is None or abs(got - target) < abs(best[1] - target):
            best = (mid, got)
        lo, hi = (mid, hi) if got < target else (lo, mid)
    return best[0]


def b80(sd, rho, n, seed):
    x = np.arange(n, dtype=float)
    lo, hi = 0.0, 1.0
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        V = ar1(R_POW, n, sd, rho, np.random.default_rng(seed)) + mid * (x - x.mean())[None, :]
        if np.mean(fit_many(V, x)[2] != 0) < 0.8:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def arma_sim(ar, ma, sig, n, reps, rng, burn=300):
    e = rng.standard_normal((reps, n + burn)) * sig
    return lfilter(np.r_[1, np.asarray(ma)], np.r_[1, -np.asarray(ar)], e, axis=1)[:, burn:]


def fitted_models(res):
    best = None
    for p in range(0, 6):
        m = ARIMA(res, order=(p, 0, 0), trend="n").fit()
        if best is None or m.aic < best[0]:
            best = (m.aic, p, m)
    _, p, mar = best
    marma = ARIMA(res, order=(1, 0, 1), trend="n").fit()
    return {f"AR({p}) by AIC": (mar.arparams if p else np.array([]), np.array([]), np.sqrt(mar.params[-1])),
            "ARMA(1,1)": (marma.arparams, marma.maparams, np.sqrt(marma.params[-1]))}


def detrended(v):
    x = np.arange(len(v), dtype=float)
    b, a = np.polyfit(x, v, 1)
    return v - (a + b * x)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "data/storm_flags.csv"
    out = sys.argv[2] if len(sys.argv) > 2 else "data"
    size1, power, sizef, calp = [], [], [], []
    for dom in ("WNP basin", "PAR"):
        fl = load_flags(src, dom)
        series = {"total": annual(fl, pd.Series(True, index=fl.index)),
                  "TS or stronger": annual(fl, any4(fl, "TS"))}
        for name, s in series.items():
            v = s.loc[1951:2023].to_numpy(float)
            res = detrended(v)
            sd, r1 = res.std(ddof=2), np.corrcoef(res[:-1], res[1:])[0, 1]
            rho = calib_rho(r1)
            print(f"{dom:10s} {name:15s} mean {v.mean():6.2f}  detrended sd {sd:.3f}  "
                  f"variance/mean {sd ** 2 / v.mean():.3f}  lag-1 {r1:+.3f}  latent rho {rho:.3f}")
            for n in NS:
                A = ar1(R_SIZE, n, sd, rho, np.random.default_rng(200 + n))
                rej = float(np.mean(fit_many(A, np.arange(n, dtype=float))[2] != 0))
                size1.append(dict(domain=dom, series=name, n=n, size_pct=100 * rej,
                                  mc_se_pct=100 * np.sqrt(rej * (1 - rej) / R_SIZE)))
            for n in NS:
                b = b80(sd, rho, n, 300 + n)
                wm = s.loc[2024 - n:2023].mean()
                power.append(dict(domain=dom, series=name, n=n, window=f"{2024 - n}-2023",
                                  window_mean=wm, b80=b, change_pct_of_mean=100 * b * (n - 1) / wm))
            for lab, (ar, ma, sig) in fitted_models(res).items():
                par = "ar=" + ",".join(f"{c:.2f}" for c in ar) + (" ma=" + ",".join(f"{c:.2f}" for c in ma) if len(ma) else "")
                for n in NS:
                    Y = arma_sim(ar, ma, sig, n, R_SIZE, np.random.default_rng(500 + n))
                    rej = float(np.mean(fit_many(Y, np.arange(n, dtype=float))[2] != 0))
                    sizef.append(dict(domain=dom, series=name, model=lab, params=par, n=n, size_pct=100 * rej,
                                      mc_se_pct=100 * np.sqrt(rej * (1 - rej) / R_SIZE)))
        cser = {"raw total": series["total"],
                "any of 3, all systems": annual(fl, requirement(fl, "any of 3", "all")),
                "JMA|JTWC, all systems": annual(fl, requirement(fl, "JMA|JTWC", "all")),
                "TS or stronger": series["TS or stronger"]}
        for name, s in cser.items():
            models = fitted_models(detrended(s.loc[1951:2023].to_numpy(float)))
            for y0 in (1951, 1961, 1977):
                vv = s.loc[y0:2023].to_numpy(float)
                xx = np.arange(len(vv), dtype=float)
                sl, band, vd = fit_many(vv[None, :], xx)
                for lab, (ar, ma, sig) in models.items():
                    Z = arma_sim(ar, ma, sig, len(vv), R_P, np.random.default_rng(900 + y0))
                    p = float(np.mean(np.abs(fit_many(Z, xx)[0]) >= abs(sl[0])))
                    calp.append(dict(domain=dom, series=name, start=y0, slope=sl[0],
                                     test_verdict={-1: "-", 0: "0", 1: "+"}[int(vd[0])],
                                     model=lab.split(" ")[0] if "AIC" in lab else lab,
                                     model_role="AIC" if "AIC" in lab else "ARMA(1,1)", p_two_sided=p))
    pd.DataFrame(size1).to_csv(os.path.join(out, "size_ar1.csv"), index=False)
    P = pd.DataFrame(power)
    P.to_csv(os.path.join(out, "power.csv"), index=False)
    F = pd.DataFrame(sizef)
    F.to_csv(os.path.join(out, "size_fitted.csv"), index=False)
    pd.DataFrame(calp).to_csv(os.path.join(out, "calibrated_p.csv"), index=False)
    pd.set_option("display.width", 200)
    print(F.pivot_table(index=["domain", "series", "model", "params"], columns="n", values="size_pct").round(1).to_string())
    print(f"written: size_ar1, power, size_fitted, calibrated_p in {out}/")


if __name__ == "__main__":
    main()
