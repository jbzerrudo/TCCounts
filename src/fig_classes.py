"""
Figure: annual main-track counts by storm class, basin and PAR, 1951-2023.

Rows are the classes of decompose.py: tropical storm or stronger, depression
only, and no agency intensity (the fourth class, below depression, is small and
is reported in the caption). Each series carries the model dating.py selects for
it: a constant, a smooth trend, or step means. Grey bands mark the periods in
which the dated steps fall. With PAGASA's annual series (from pagasa.py) the PAR
panels also show PAGASA's own count, classified by its 6-hourly categories
inside the PAR, 1960-2023. Without it the figure is drawn from IBTrACS alone.

Usage
    python src/fig_classes.py data/class_counts.csv figures/ [pagasa_private/pagasa_annual.csv]
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

from dating import Poisson, date_series

MM = 1 / 25.4
INK, PAG = "#333333", "#008300"           # IBTrACS; PAGASA (palette slot 6, clear of the agency colours)
BAND, AXIS, TEXT2 = "#ECECEA", "#898781", "#52514E"
PERIODS = [(1967, 1972), (1990, 1993), (2003, 2004)]
ROWS = [("ts_or_stronger", "pagasa_ts", "tropical storm or stronger"),
        ("depression_only", "pagasa_td", "depression only"),
        ("no_intensity", None, "graded by no agency")]
DOMAINS = [("WNP basin", "Basin"), ("PAR", "PAR")]

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": TEXT2, "ytick.labelcolor": TEXT2, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})


def fitted(v, t):
    """The model dating.py selects for this series, evaluated at every year."""
    m = Poisson(np.asarray(v, float))
    r = date_series(m, t)
    if r["selected"] == "smooth":
        return m.smooth(t)[0], r
    br = tuple(int(np.flatnonzero(t == int(y))[0]) for y in str(r["breaks"]).split()) if r["breaks"] else ()
    return m.step(br), r


def main():
    cc, out = sys.argv[1], sys.argv[2]
    pg = pd.read_csv(sys.argv[3], index_col=0) if len(sys.argv) > 3 else None
    os.makedirs(out, exist_ok=True)
    c = pd.read_csv(cc)
    t = np.arange(1951, 2024)
    fig, axes = plt.subplots(3, 2, figsize=(174 * MM, 150 * MM), sharex=True)
    letters = iter("abcdef")
    report = []
    for i, (cls, pcol, label) in enumerate(ROWS):
        for j, (dom, short) in enumerate(DOMAINS):
            ax = axes[i, j]
            v = c[c.domain == dom].set_index("SEASON")[cls].loc[1951:2023].to_numpy(float)
            for a, b in PERIODS:
                ax.axvspan(a - 0.5, b + 0.5, color=BAND, lw=0, zorder=0)
            f, r = fitted(v, t)
            ax.plot(t, v, color=INK, lw=0.6, alpha=0.45, zorder=2)
            ax.step(t, f, where="mid", color=INK, lw=1.6, zorder=3)
            report.append((short, cls, r["selected"], r["breaks"], r["segments"]))
            top = v.max()
            if pg is not None and pcol is not None and dom == "PAR":
                tp = pg.index.to_numpy()
                pv = pg[pcol].to_numpy(float)
                fp, rp = fitted(pv, tp)
                ax.plot(tp, pv, color=PAG, lw=0.6, alpha=0.6, zorder=2)
                ax.step(tp, fp, where="mid", color=PAG, lw=1.6, zorder=3)
                report.append((short, pcol, rp["selected"], rp["breaks"], rp["segments"]))
                top = max(top, pv.max())
            ax.set_ylim(0, top * 1.12)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True, steps=[1, 2, 5, 10]))
            ax.set_xlim(1950, 2024)
            ax.set_xticks(range(1950, 2030, 10))
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.set_title(f"({next(letters)}) {short}, {label}", loc="left", fontsize=9, color="#0B0B0B", pad=4)
            if j == 0:
                ax.set_ylabel("Storms per year", color=TEXT2)
            if i == 2:
                ax.set_xlabel("Year", color=TEXT2)
    handles = [Line2D([], [], color=INK, lw=0.6, alpha=0.45, label="IBTrACS, annual count"),
               Line2D([], [], color=INK, lw=1.6, label="IBTrACS, selected model"),
               Line2D([], [], color=PAG, lw=0.6, alpha=0.6, label="PAGASA, annual count"),
               Line2D([], [], color=PAG, lw=1.6, label="PAGASA, selected model"),
               Patch(color=BAND, lw=0, label="Periods with dated steps")]
    if pg is None:
        handles = handles[:2] + handles[4:]
    fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, 1.0), handlelength=2.2, columnspacing=1.6, labelcolor=TEXT2)
    fig.tight_layout(rect=(0, 0, 1, 0.93), h_pad=1.0, w_pad=1.5)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"Fig3_classes.{ext}"), dpi=300)
    for row in report:
        print(" | ".join(str(x) for x in row))
    print(f"written: {out}/Fig3_classes.pdf and .png")


if __name__ == "__main__":
    main()
