"""Regenerate data/start_year_sweep.csv from the annual series."""
import sys
import numpy as np, pandas as pd, scipy.stats as st

SER = sys.argv[1] if len(sys.argv) > 1 else "data/par_annual_series.csv"
DST = sys.argv[2] if len(sys.argv) > 2 else "data/start_year_sweep.csv"
S = pd.read_csv(SER).set_index("SEASON")
rng = np.random.default_rng(41)


def sweep_one(ser, L=10, nb=6000):
    v = np.asarray(ser, float); x = np.asarray(ser.index, float)
    xc = x - x.mean(); den = (xc ** 2).sum(); N = len(v)
    sl = (v @ xc) / den
    starts = rng.integers(0, N, size=(nb, int(np.ceil(N / L))))
    idx = (starts[:, :, None] + np.arange(L)[None, None, :]).reshape(nb, -1)[:, :N] % N
    lo, hi = np.percentile((v[idx] @ xc) / den, [2.5, 97.5])
    return sl, lo, hi, ("+" if sl > hi else "-" if sl < lo else "0")


rows = []
for y in range(1884, 2001, 2):
    sl, lo, hi, sig = sweep_one(S.loc[y:2023, "total"])
    if y >= 1951:
        cs, clo, chi, cg = sweep_one(S.loc[y:2023, "any_int"])
    else:
        cs = clo = chi = np.nan; cg = np.nan
    rows.append((y, sl, lo, hi, sig, cs, clo, chi, cg))
pd.DataFrame(rows, columns=["start", "slope", "lo", "hi", "sig",
                            "cslope", "clo", "chi", "csig"]).to_csv(DST, index=False)
print(f"written: {DST}")