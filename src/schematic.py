"""
Fig 2: the enumeration. What each of the three choices does to the count, and
how 2,400 specifications arise.

Purely structural: no data are read, so this figure is deterministic.

Usage
    python src/schematic.py figures/
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

INK, MUTE, ACCENT = "#1A1A1A", "#6E6E6E", "#1C635F"
FILL_STAGE, FILL_CHOICE = "#F2F2F0", "#FFFFFF"
W, XC = 7.2, 0.40           # figure width; centre of the flow column, in axes units

STAGES = [
    ("src",    "IBTrACS v04r01, western North Pacific",
               "4,076 main-track storms, 1884 to 2023", None, None),
    ("choice", "domain",
               "WNP basin   |   PAR hexagon", "2", "2"),
    ("choice", "counting rule",
               "all SIDs  |  main track  |  $\\geq$2 fixes  |  tropical nature", "4", "8"),
    ("choice", "grading requirement",
               "raw total  |  JTWC  |  JMA  |  JMA+JTWC  |  JMA+JTWC+CMA  |  all four", "6", "48"),
    ("proc",   "annual count series",
               "one count per season, 1884 to 2023", None, None),
    ("choice", "start year",
               "1951 to 2000, every window ending 2023", "50", "2,400"),
    ("proc",   "ordinary least squares slope",
               "standard error corrected for serial dependence", None, None),
    ("out",    "verdict",
               "significant decline  |  no detectable trend  |  significant increase", None, None),
]


def main():
    dst = (sys.argv[1] if len(sys.argv) > 1 else "figures/").rstrip("/")
    n = len(STAGES)
    fig, ax = plt.subplots(figsize=(W, 0.86 * n + 0.5))
    ax.set_xlim(0, 1); ax.set_ylim(0, n); ax.axis("off")

    bw, bh = 0.66, 0.52
    for i, (kind, title, sub, mult, run) in enumerate(STAGES):
        y = n - 0.65 - i
        choice = kind == "choice"
        ax.add_patch(FancyBboxPatch(
            (XC - bw / 2, y - bh / 2), bw, bh,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            linewidth=1.5 if choice else 0.9,
            edgecolor=ACCENT if choice else MUTE,
            facecolor=FILL_CHOICE if choice else FILL_STAGE, zorder=3))
        ax.text(XC, y + 0.115, title, ha="center", va="center", zorder=4,
                fontsize=10 if choice else 9.4, color=ACCENT if choice else INK,
                fontweight="bold" if choice else "normal")
        ax.text(XC, y - 0.125, sub, ha="center", va="center", zorder=4,
                fontsize=7.6, color=MUTE if not choice else INK)

        if i < n - 1:
            ax.add_patch(FancyArrowPatch(
                (XC, y - bh / 2), (XC, y - 1 + bh / 2),
                arrowstyle="-|>", mutation_scale=11, linewidth=1.0,
                color=MUTE, shrinkA=0, shrinkB=0, zorder=2))

        if mult:
            ax.text(XC + bw / 2 + 0.035, y + 0.10, r"$\times$" + mult,
                    ha="left", va="center", fontsize=9.6, color=ACCENT, fontweight="bold")
            ax.text(XC + bw / 2 + 0.035, y - 0.115, run + " so far",
                    ha="left", va="center", fontsize=7.4, color=MUTE)

    fig.tight_layout(pad=0.3)
    for ext in ("pdf", "png"):
        fig.savefig(f"{dst}/Fig2_design.{ext}", dpi=400, bbox_inches="tight")
    print(f"written: {dst}/Fig2_design.pdf")


if __name__ == "__main__":
    main()
