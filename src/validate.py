"""
Establish, by simulation, the size and power of the significance test used in
multiverse.py. These are the numbers quoted in Section 2 of the manuscript.

Usage
    python src/validate.py            # ~3 minutes
"""
import numpy as np
from scipy import stats

REPS_SIZE, REPS_POWER = 3000, 600


def verdict(v, x):
    """The test from multiverse.py: OLS slope against a serial-dependence
    corrected standard error."""
    n = len(v)
    xc = x - x.mean()
    den = (xc ** 2).sum()
    slope = (v @ xc) / den
    resid = v - np.poly1d(np.polyfit(x, v, 1))(x)
    r1 = min(max(np.corrcoef(resid[:-1], resid[1:])[0, 1], 0.0), 0.99)
    n_eff = max(n * (1 - r1) / (1 + r1), 4.0)
    se = np.sqrt((resid @ resid) / (n - 2) / den) * np.sqrt(n / n_eff)
    return abs(slope / se) > stats.t.ppf(0.975, max(n_eff - 2, 1))


def ar1_counts(n, mu, rho, rng):
    e = rng.standard_normal(n)
    z = np.empty(n); z[0] = e[0]
    for i in range(1, n):
        z[i] = rho * z[i - 1] + np.sqrt(1 - rho ** 2) * e[i]
    return rng.poisson(np.maximum(mu + np.sqrt(mu) * z, 0.1)).astype(float)


def size():
    print("\nEmpirical size, nominal 5 percent, %d realisations" % REPS_SIZE)
    print(f"{'lag-1 autocorrelation':<24}{'n=73':>8}{'n=33':>8}")
    for rho in (0.0, 0.15, 0.33, 0.50):
        row = []
        for n in (73, 33):
            x = np.arange(n, dtype=float)
            hits = sum(verdict(ar1_counts(n, 29.1, rho,
                                          np.random.default_rng(510000 + k * 17 + int(rho * 1000) + n)), x)
                       for k in range(REPS_SIZE))
            row.append(100 * hits / REPS_SIZE)
        print(f"{rho:<24.2f}" + "".join(f"{r:>7.1f}%" for r in row))


def power():
    print("\nPower at n=73, %d realisations" % REPS_POWER)
    print(f"{'mean':>6}{'trend':>9}{'total change':>14}{'% of mean':>11}{'power':>8}")
    for mu in (19.7, 29.1):
        for beta in (0.04, 0.06, 0.07, 0.08, 0.10):
            n = 73
            x = np.arange(n, dtype=float)
            hits = 0
            for k in range(REPS_POWER):
                rng = np.random.default_rng(620000 + k * 13 + int(mu * 10) + int(beta * 1000))
                lam = np.maximum(mu + beta * (x - x.mean()), 0.1)
                hits += verdict(rng.poisson(lam).astype(float), x)
            tot = beta * (n - 1)
            print(f"{mu:>6.1f}{beta:>+9.2f}{tot:>+14.1f}{100*tot/mu:>10.0f}%{100*hits/REPS_POWER:>7.1f}%")


if __name__ == "__main__":
    size()
    power()
