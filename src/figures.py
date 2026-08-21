"""
Regenerate both manuscript figures from the archived data.

    Fig 1  counts and agency coverage, basin and PAR
    Fig 2  the specification curve over all 2,400 fits

Usage
    python src/figures.py par_clipped.csv ibtracs.WP.list.v04r01.csv \
                          data/multiverse.csv figures/
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FIRST, LAST = 1884, 2023
YEARS = np.arange(FIRST, LAST + 1)
SIG, NUL = "#D55E00", "#9A9A9A"          # status pair: significant / not
AGENCIES = ["raw total", "all four", "JMA+JTWC+CMA", "JMA", "JMA+JTWC", "JTWC"]
RULES = ["all SIDs", "main track", ">=2 fixes", "tropical nature"]


def num(d, c):
    return pd.to_numeric(d[c], errors="coerce").fillna(0)


def coverage(d):
    """Annual storm count, and the count each agency assigned an intensity to."""
    d = d[(d.SEASON >= FIRST) & (d.SEASON <= LAST)]
    m = d[d.TRACK_TYPE == "main"]
    n = lambda sub: sub.groupby("SEASON").SID.nunique().reindex(YEARS).fillna(0).astype(int)
    gcol = "TOKYO_GRADE" if "TOKYO_GRADE" in m.columns else "TOK_GRADE"
    cov = {"HKO": n(m[num(m, "HKO_WIND") > 0]),
           "JMA": n(m[num(m, gcol).isin([2, 3, 4, 5, 9])]),
           "CMA": n(m[num(m, "CMA_WIND") > 0]),
           "JTWC": n(m[num(m, "USA_WIND") > 0])}
    return n(m).values, {k: v.values for k, v in cov.items()}


def fig1(dom, out):
    ramp = plt.get_cmap("Blues")(np.linspace(0.40, 0.95, 4))
    fig, ax = plt.subplots(2, 2, figsize=(10.4, 5.7))
    for j, name in enumerate(["WNP basin", "PAR"]):
        tot, cov = dom[name]
        a = ax[0, j]
        a.fill_between(YEARS, tot, color="#DFDFDF", lw=0)
        a.plot(YEARS, tot, color="#333333", lw=.85)
        a.set_xlim(FIRST, LAST); a.set_ylim(0, 58 if j == 0 else 40)
        a.set_title(f"({'ab'[j]})  {name}, all track entries", loc="left", fontsize=9.5)
        a.set_xlabel("Year", fontsize=9)
        if j == 0:
            a.set_ylabel("Storms per year")
        b = ax[1, j]
        for c, k in zip(ramp, ["HKO", "JMA", "CMA", "JTWC"]):
            share = 100 * np.where(tot > 0, cov[k] / np.maximum(tot, 1), np.nan)
            sm = np.convolve(share, np.ones(9) / 9, mode="same")
            b.plot(YEARS[4:-4], sm[4:-4], color=c, lw=1.9, label=k)
        b.set_ylim(0, 105); b.set_xlim(1944, LAST)
        b.set_title(f"({'cd'[j]})  share of those storms an agency graded",
                    loc="left", fontsize=9.5)
        b.set_xlabel("Year", fontsize=9)
        if j == 0:
            b.set_ylabel("Percent")
            b.legend(fontsize=7.8, frameon=False, loc="lower right", ncol=2,
                     handlelength=1.4, columnspacing=1.1)
        for p in (a, b):
            for s in ("top", "right"):
                p.spines[s].set_visible(False)
            p.tick_params(labelsize=8)
    plt.tight_layout(pad=.6, w_pad=1.9, h_pad=1.4)
    for ext in ("pdf", "png"):
        plt.savefig(f"{out}Fig1.{ext}", dpi=190, bbox_inches="tight")
    plt.close()


def fig2(R, out):
    """Fitted slope against start year, one panel per agency set, both domains."""
    R = R[R.rule == "main track"]
    colour = {"WNP basin": "#0072B2", "PAR": "#D55E00"}
    fig, ax = plt.subplots(2, 3, figsize=(10.6, 5.8), sharex=True, sharey=True)
    for i, a in enumerate(AGENCIES):
        p = ax[i // 3, i % 3]
        for dom, c in colour.items():
            g = R[(R.agency == a) & (R.domain == dom)].sort_values("start")
            p.plot(g.start, g.lo, color=c, lw=.7, ls=(0, (3, 2)), alpha=.75)
            p.plot(g.start, g.hi, color=c, lw=.7, ls=(0, (3, 2)), alpha=.75)
            p.plot(g.start, g.slope, color=c, lw=1.9, label=dom, zorder=4)
            sig = g[g.verdict == "-"]
            p.plot(sig.start, sig.slope, "o", ms=4.0, color=c, mec="white", mew=.7, zorder=5)
        p.axhline(0, color="#666666", lw=.8)
        p.set_title(a, loc="left", fontsize=9.5)
        p.set_ylim(-0.40, 0.24); p.set_xlim(1951, 2000)
        for s in ("top", "right"):
            p.spines[s].set_visible(False)
        p.tick_params(labelsize=8)
        if i % 3 == 0:
            p.set_ylabel("OLS slope (storms yr$^{-1}$)", fontsize=9)
        if i // 3 == 1:
            p.set_xlabel("First year of the fitted window", fontsize=9)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend(h, l, fontsize=8.6, frameon=False, ncol=2, loc="upper left",
               bbox_to_anchor=(0.075, 1.055), handlelength=1.6, columnspacing=1.6)
    fig.text(0.40, 1.032,
             "dashed: no-trend band          filled circles: significant decline",
             fontsize=8.2, color="#444444")
    plt.tight_layout(pad=.6, w_pad=1.3, h_pad=1.2)
    for ext in ("pdf", "png"):
        plt.savefig(f"{out}Fig2.{ext}", dpi=185, bbox_inches="tight")
    plt.close()


def main():
    par = sys.argv[1] if len(sys.argv) > 1 else "par_clipped.csv"
    wnp = sys.argv[2] if len(sys.argv) > 2 else "ibtracs.WP.list.v04r01.csv"
    mvf = sys.argv[3] if len(sys.argv) > 3 else "data/multiverse.csv"
    out = sys.argv[4] if len(sys.argv) > 4 else "figures/"
    dom = {"WNP basin": coverage(pd.read_csv(wnp, low_memory=False, skiprows=[1])),
           "PAR": coverage(pd.read_csv(par, low_memory=False))}
    fig1(dom, out)
    fig2(pd.read_csv(mvf), out)
    print(f"written: {out}Fig1.pdf, {out}Fig2.pdf")


if __name__ == "__main__":
    main()
