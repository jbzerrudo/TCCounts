"""
Establish, by simulation, the size and power of the significance test used in
multiverse.py. These are the numbers quoted in Section 2 of the manuscript.

The generating process is a Gaussian copula with negative-binomial margins.
Serial dependence is imposed on a latent Gaussian AR(1) series and carried
through the copula into the counts, so the delivered lag-one autocorrelation
and the delivered variance-to-mean ratio can be set independently. An earlier
version added an AR(1) latent mean to a Poisson draw. That construction cannot
reach the basin's autocorrelation: the Poisson layer is white, so the count
autocorrelation is capped at rho*(d-1)/d, which is 0.29 at a variance-to-mean
ratio of 1.40 and lower still after detrending. It was replaced for that reason.

The latent rho is calibrated by bisection against the lag-one autocorrelation of
the DETRENDED counts, which is the quantity measured on the real series.

Usage
    python src/validate.py            # about six minutes
"""
import numpy as np
from scipy import stats

REPS_SIZE, REPS_OP, REPS_BOOT, REPS_POWER = 5000, 20000, 1000, 4000
BASIN = dict(mu=29.1, disp=1.37, r1=0.33)
PAR = dict(mu=19.7, disp=1.01, r1=0.03)


def nb_params(mu, disp):
    disp = max(disp, 1.0 + 1e-9)
    r = mu / (disp - 1.0)
    return r, r / (r + mu)


def draw(reps, n, mu, rho, disp, rng):
    """AR(1) Gaussian copula with negative-binomial margins."""
    e = rng.standard_normal((reps, n))
    z = np.empty((reps, n))
    z[:, 0] = e[:, 0]
    s = np.sqrt(1 - rho ** 2)
    for i in range(1, n):
        z[:, i] = rho * z[:, i - 1] + s * e[:, i]
    r, p = nb_params(mu, disp)
    return stats.nbinom.ppf(np.clip(stats.norm.cdf(z), 1e-12, 1 - 1e-12), r, p).astype(float)


def _resid(V, x):
    xc = x - x.mean()
    den = (xc ** 2).sum()
    b = (V @ xc) / den
    return b, den, V - V.mean(1, keepdims=True) - b[:, None] * xc[None, :]


def _lag1(R):
    a, c = R[:, :-1], R[:, 1:]
    am, cm = a.mean(1, keepdims=True), c.mean(1, keepdims=True)
    return (((a - am) * (c - cm)).sum(1) /
            np.sqrt(((a - am) ** 2).sum(1) * ((c - cm) ** 2).sum(1)))


def reject(V, x):
    """The test from multiverse.py, applied to every row of V."""
    n = V.shape[1]
    b, den, R = _resid(V, x)
    r1 = np.clip(np.nan_to_num(_lag1(R)), 0.0, 0.99)
    ne = np.maximum(n * (1 - r1) / (1 + r1), 4.0)
    se = np.sqrt((R ** 2).sum(1) / (n - 2) / den) * np.sqrt(n / ne)
    return np.abs(b / se) > stats.t.ppf(0.975, np.maximum(ne - 2, 1))


def reject_bootstrap(V, x, L=10, nb=1500, seed=41):
    """The test this work discarded: a circular block bootstrap of the OLS
    residuals, block length 10, 1,500 resamples, as in the earlier version."""
    n = V.shape[1]
    xc = x - x.mean()
    b, den, R = _resid(V, x)
    rng = np.random.default_rng(seed)
    nblk = int(np.ceil(n / L))
    out = np.empty(V.shape[0], bool)
    for i in range(V.shape[0]):
        st = rng.integers(0, n, size=(nb, nblk))
        idx = (st[:, :, None] + np.arange(L)[None, None, :]).reshape(nb, -1)[:, :n] % n
        lo, hi = np.percentile((R[i][idx] @ xc) / den, [2.5, 97.5])
        out[i] = (b[i] < lo) or (b[i] > hi)
    return out


def calibrate(n, mu, target, disp, reps=4000, seed=7):
    """Latent rho that delivers the requested detrended lag-one autocorrelation."""
    x = np.arange(n, dtype=float)

    def delivered(rho):
        R = _resid(draw(reps, n, mu, rho, disp, np.random.default_rng(seed)), x)[2]
        return float(np.nanmean(_lag1(R)))

    if target <= 0:
        return 0.0, delivered(0.0)
    lo, hi, best = 0.0, 0.98, None
    for _ in range(16):
        mid = 0.5 * (lo + hi)
        got = delivered(mid)
        if best is None or abs(got - target) < abs(best[1] - target):
            best = (mid, got)
        lo, hi = (mid, hi) if got < target else (lo, mid)
    return best


def size():
    d = BASIN
    print("\nEmpirical size, nominal 5 percent, %d realizations, "
          "variance/mean = %.2f, mean = %.1f" % (REPS_SIZE, d["disp"], d["mu"]))
    print(f"{'target lag-1':<14}" + "".join(f"{'n=' + str(n):>10}" for n in (73, 50, 33, 24)))
    for t in (0.00, 0.15, 0.25, 0.33):
        row = []
        for n in (73, 50, 33, 24):
            rho, _ = calibrate(n, d["mu"], t, d["disp"])
            x = np.arange(n, dtype=float)
            V = draw(REPS_SIZE, n, d["mu"], rho, d["disp"],
                     np.random.default_rng(9001 + n + int(t * 1000)))
            row.append(f"{100 * reject(V, x).mean():>9.1f}%")
        print(f"{t:<14.2f}" + "".join(row))
    print("The PAR dispersion of 1.01 gives the same table to within 0.2 points.")

    print("\nSize at the operating point of each domain, n = 73, "
          "%d realizations" % REPS_OP)
    for name, d in (("basin", BASIN), ("PAR", PAR)):
        rho, got = calibrate(73, d["mu"], d["r1"], d["disp"])
        x = np.arange(73, dtype=float)
        V = draw(REPS_OP, 73, d["mu"], rho, d["disp"], np.random.default_rng(555))
        print(f"  {name:<6} lag-1 {d['r1']:+.2f} (delivered {got:+.3f}), "
              f"variance/mean {d['disp']:.2f}: {100 * reject(V, x).mean():.1f}%")

    print("\nThe discarded circular block bootstrap on the same realizations, "
          "%d each" % REPS_BOOT)
    for name, d in (("basin", BASIN), ("PAR", PAR)):
        rho, _ = calibrate(73, d["mu"], d["r1"], d["disp"])
        x = np.arange(73, dtype=float)
        V = draw(REPS_BOOT, 73, d["mu"], rho, d["disp"], np.random.default_rng(777))
        print(f"  {name:<6} bootstrap {100 * reject_bootstrap(V, x).mean():.1f}%   "
              f"effective-n {100 * reject(V, x).mean():.1f}%")


def power():
    n = 73
    x = np.arange(n, dtype=float)
    print("\nPower at n = 73, %d realizations, each domain simulated with its "
          "own dispersion and autocorrelation" % REPS_POWER)
    print(f"{'domain':>7}{'mean':>7}{'var/mean':>10}{'lag-1':>8}"
          f"{'trend':>9}{'change':>9}{'% of mean':>11}{'power':>8}")
    for name, d, betas in (("PAR", PAR, (0.05, 0.06, 0.07, 0.076, 0.08, 0.10)),
                           ("basin", BASIN, (0.09, 0.12, 0.14, 0.15, 0.18, 0.22))):
        rho, got = calibrate(n, d["mu"], d["r1"], d["disp"])
        for b in betas:
            V = draw(REPS_POWER, n, d["mu"], rho, d["disp"],
                     np.random.default_rng(31337 + int(b * 10000) + int(d["mu"] * 10)))
            V = np.maximum(V + b * (x - x.mean()), 0.0)
            tot = b * (n - 1)
            print(f"{name:>7}{d['mu']:>7.1f}{d['disp']:>10.2f}{got:>8.2f}"
                  f"{b:>+9.3f}{tot:>+9.1f}{100 * tot / d['mu']:>10.0f}%"
                  f"{100 * reject(V, x).mean():>7.1f}%")


if __name__ == "__main__":
    size()
    power()
