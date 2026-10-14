"""
Figure: the enumeration. The four choices that enter a fitted count trend in
this study and how the specifications arise. Purely structural: no data are
read, so the figure is deterministic. Adapted from schematic.py (version 3),
with the intensity threshold added as a fourth choice.

Usage
    python src/fig_design.py figures/
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

MM = 1 / 25.4
INK, MUTE, ACCENT = "#1A1A1A", "#6E6E6E", "#1C635F"
FILL_STAGE, FILL_CHOICE = "#F2F2F0", "#FFFFFF"
XC = 0.40

STAGES = [
    ("src", "IBTrACS v04r01, western North Pacific", "4,076 main-track storms, 1884 to 2023", None, None),
    ("choice", "domain", "WNP basin  |  PAR hexagon", "2", "2"),
    ("choice", "counting rule", "all SIDs  |  main track  |  $\\geq$2 fixes  |  tropical nature", "4", "8"),
    ("choice", "grading requirement", "JMA  |  JTWC  |  CMA  |  JMA or JTWC  |  any of the three  |  all three", "6", "48"),
    ("choice", "intensity threshold", "all systems  |  depression or stronger  |  tropical storm or stronger", "3", "144"),
    ("proc", "annual count series", "one count per season, 1884 to 2023", None, None),
    ("choice", "start year", "1951 to 2000, every window ending in 2023", "50", "7,200"),
    ("proc", "ordinary least squares slope", "standard error corrected for serial dependence", None, None),
    ("out", "verdict", "significant decline  |  no detectable trend  |  significant increase", None, None),
]


def main():
    dst = sys.argv[1] if len(sys.argv) > 1 else "figures"
    os.makedirs(dst, exist_ok=True)
    plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
    n = len(STAGES)
    fig, ax = plt.subplots(figsize=(174 * MM, 0.66 * n + 0.3))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")
    bw, bh = 0.70, 0.56
    for i, (kind, title, sub, mult, run) in enumerate(STAGES):
        y = n - 0.62 - i
        choice = kind == "choice"
        ax.add_patch(FancyBboxPatch((XC - bw / 2, y - bh / 2), bw, bh,
                                    boxstyle="round,pad=0.012,rounding_size=0.02",
                                    linewidth=1.4 if choice else 0.8,
                                    edgecolor=ACCENT if choice else MUTE,
                                    facecolor=FILL_CHOICE if choice else FILL_STAGE, zorder=3))
        ax.text(XC, y + 0.12, title, ha="center", va="center", zorder=4, fontsize=9.5 if choice else 9,
                color=ACCENT if choice else INK, fontweight="bold" if choice else "normal")
        ax.text(XC, y - 0.13, sub, ha="center", va="center", zorder=4, fontsize=8,
                color=INK if choice else MUTE)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((XC, y - bh / 2), (XC, y - 1 + bh / 2), arrowstyle="-|>",
                                         mutation_scale=10, linewidth=0.9, color=MUTE, shrinkA=0,
                                         shrinkB=0, zorder=2))
        if mult:
            ax.text(XC + bw / 2 + 0.03, y + 0.11, r"$\times$" + mult, ha="left", va="center",
                    fontsize=9.5, color=ACCENT, fontweight="bold")
            ax.text(XC + bw / 2 + 0.03, y - 0.13, run + " so far", ha="left", va="center",
                    fontsize=8, color=MUTE)
    fig.tight_layout(pad=0.2)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(dst, f"ESM_FigS1_design.{ext}"), dpi=300)
    print(f"written: {dst}/ESM_FigS1_design.pdf")


if __name__ == "__main__":
    main()
