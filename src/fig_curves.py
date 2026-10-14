"""
Figure: fitted slope against the first year of the window, for the six grading
requirements at two intensity thresholds, main-track counting rule, basin and
PAR. Every window ends in 2023. Dashed lines are the no-trend band of the test
(multiverse.fit); filled circles mark significant declines and triangles
significant increases.

Usage
    python src/fig_curves.py data/grid_threshold.csv figures/
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

MM = 1 / 25.4
DOM = {"WNP basin": ("#0072B2", "Basin"), "PAR": ("#D55E00", "PAR")}
REQ = [("JMA", "JMA"), ("JTWC", "JTWC"), ("CMA", "CMA"), ("JMA|JTWC", "JMA or JTWC"),
       ("any of 3", "Any of the three"), ("all 3", "All three")]
THR = [("all", "all systems"), ("TS", "tropical storm or stronger")]
AXIS, TEXT2 = "#898781", "#52514E"

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": TEXT2, "ytick.labelcolor": TEXT2, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
})


def main():
    src, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    G = pd.read_csv(src)
    G = G[G.rule == "main track"]
    fig, axes = plt.subplots(len(REQ), 2, figsize=(174 * MM, 222 * MM), sharex=True, sharey=True)
    letters = iter("abcdefghijkl")
    for i, (req, rname) in enumerate(REQ):
        for j, (thr, tname) in enumerate(THR):
            ax = axes[i, j]
            ax.axhline(0, color="#B5B4AD", lw=0.6, zorder=1)
            for dom, (col, _) in DOM.items():
                g = G[(G.domain == dom) & (G.requirement == req) & (G.threshold == thr)].sort_values("start")
                ax.plot(g.start, g.lo, color=col, lw=0.6, ls=(0, (3, 2)), alpha=0.8, zorder=2)
                ax.plot(g.start, g.hi, color=col, lw=0.6, ls=(0, (3, 2)), alpha=0.8, zorder=2)
                ax.plot(g.start, g.slope, color=col, lw=1.3, zorder=3)
                d = g[g.verdict == "-"]
                ax.plot(d.start, d.slope, "o", ms=3.2, color=col, mec="white", mew=0.5, zorder=4)
                u = g[g.verdict == "+"]
                ax.plot(u.start, u.slope, "^", ms=4.5, color=col, mec="white", mew=0.5, zorder=5)
            ax.set_title(f"({next(letters)}) {rname}, {tname}", loc="left", fontsize=8.5, color="#0B0B0B", pad=3)
            ax.set_xlim(1949.5, 2001.5)
            ax.set_ylim(-0.36, 0.36)
            ax.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            if i == len(REQ) - 1:
                ax.set_xlabel("First year of the fitted window", color=TEXT2)
    handles = [Line2D([], [], color=DOM["WNP basin"][0], lw=1.3, label="Basin"),
               Line2D([], [], color=DOM["PAR"][0], lw=1.3, label="PAR"),
               Line2D([], [], color=AXIS, lw=0.6, ls=(0, (3, 2)), label="No-trend band"),
               Line2D([], [], color=AXIS, lw=0, marker="o", ms=3.2, label="Significant decline"),
               Line2D([], [], color=AXIS, lw=0, marker="^", ms=4.5, label="Significant increase")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=8,
               bbox_to_anchor=(0.5, 1.0), handlelength=2.0, columnspacing=1.4, labelcolor=TEXT2)
    fig.supylabel("Fitted slope (storms per year)", color=TEXT2, fontsize=8, x=0.01)
    fig.tight_layout(rect=(0.02, 0, 1, 0.975), h_pad=0.6, w_pad=1.2)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"Fig2_curves.{ext}"), dpi=300)
    print(f"written: {out}/Fig2_curves.pdf and .png")


if __name__ == "__main__":
    main()
