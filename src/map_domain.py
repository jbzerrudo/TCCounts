"""
Fig 1: the study area. The western North Pacific archive as track density, with
the PAGASA Philippine Area of Responsibility drawn as the hexagon it is.

Density is the number of DISTINCT main-track storms whose track passes through
each one-degree cell, not the number of best-track fixes. Fix frequency changed
across the record, so a fix-count density would carry the same bookkeeping
artifact this paper is about.

Coastline comes from data/coastline_wnp.npz, extracted once from the GSHHS
shoreline distributed with basemap-data, so this script needs only matplotlib
and numpy and runs with no network access.

Usage
    python src/map_domain.py ibtracs.WP.list.v04r01.csv figures/
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

LON0, LON1, LAT0, LAT1 = 100, 180, 0, 55
PAR = [(115, 5), (115, 15), (120, 21), (120, 25), (135, 25), (135, 5)]
PAR_RED = "#E31A1C"                            # PAR boundary
CMAP = LinearSegmentedColormap.from_list(
    "wnp", ["#FFFFFF", "#DCECE8", "#A9D3CC", "#6BB3AA", "#3B8C85", "#1C635F", "#0B3D3B"])


def coastline(path):
    z = np.load(path, allow_pickle=True)
    pts, lens = z["pts"], z["lens"]
    out, i = [], 0
    for n in lens:
        out.append(pts[i:i + n]); i += n
    return out


def density(src):
    d = pd.read_csv(src, low_memory=False, skiprows=[1],
                    usecols=["SID", "SEASON", "LAT", "LON", "TRACK_TYPE"])
    d = d[(d.SEASON >= 1884) & (d.SEASON <= 2023) & (d.TRACK_TYPE == "main")]
    lat = pd.to_numeric(d.LAT, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(d.LON, errors="coerce").to_numpy(float) % 360
    sid = d.SID.to_numpy()
    good = np.isfinite(lat) & np.isfinite(lon)
    lat, lon, sid = lat[good], lon[good], sid[good]
    nx, ny = LON1 - LON0, LAT1 - LAT0
    ix = np.floor(lon - LON0).astype(int)
    iy = np.floor(lat - LAT0).astype(int)
    m = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
    pairs = pd.DataFrame({"c": ix[m] * 1000 + iy[m], "s": sid[m]}).drop_duplicates()
    cnt = pairs.groupby("c").size()
    H = np.zeros((ny, nx))
    H[(cnt.index.values % 1000).astype(int), (cnt.index.values // 1000).astype(int)] = cnt.values
    return np.ma.masked_where(H == 0, H), d.SID.nunique(), len(lat), float(m.mean())


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "ibtracs.WP.list.v04r01.csv"
    dst = (sys.argv[2] if len(sys.argv) > 2 else "figures/").rstrip("/")
    coast = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         os.pardir, "data", "coastline_wnp.npz")
    H, n_storm, n_fix, frac = density(src)
    print(f"{n_storm} main-track storms, {n_fix} fixes, "
          f"{100 * frac:.1f}% of fixes inside the frame, "
          f"densest cell {int(H.max())} storms")

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    im = ax.pcolormesh(np.arange(LON0, LON1 + 1), np.arange(LAT0, LAT1 + 1), H,
                       cmap=CMAP, vmin=0, vmax=H.max(), zorder=2, rasterized=True)
    for seg in coastline(coast):
        ax.plot(seg[:, 0], seg[:, 1], color="#3A3A3A", lw=.55, zorder=3,
                solid_joinstyle="round", solid_capstyle="round")
    for v in range(LON0 + 10, LON1, 10):
        ax.axvline(v, color="#FFFFFF", lw=.35, alpha=.45, zorder=4)
    for v in range(LAT0 + 10, LAT1, 10):
        ax.axhline(v, color="#FFFFFF", lw=.35, alpha=.45, zorder=4)

    px = [p[0] for p in PAR] + [PAR[0][0]]
    py = [p[1] for p in PAR] + [PAR[0][1]]
    ax.plot(px, py, color=PAR_RED, lw=1.1, zorder=6, solid_joinstyle="miter")

    ax.text(139.4, 27.4, "PAR", fontsize=12.5, fontweight="bold", color="#111111",
            ha="center", va="center", zorder=8,
            bbox=dict(boxstyle="round,pad=0.30", fc="#FFFFFF", ec=PAR_RED, lw=1.0))

    ax.set_xlim(LON0, LON1); ax.set_ylim(LAT0, LAT1)
    ax.set_aspect(1.0)
    ax.set_xticks(range(LON0, LON1 + 1, 10))
    ax.set_yticks(range(LAT0, LAT1 + 1, 10))
    ax.set_xticklabels([f"{v}$^{{\\circ}}$E" for v in range(LON0, LON1 + 1, 10)], fontsize=8)
    ax.set_yticklabels([f"{v}$^{{\\circ}}$N" for v in range(LAT0, LAT1 + 1, 10)], fontsize=8)
    ax.tick_params(length=3, width=.7, color="#444444")
    for s in ax.spines.values():
        s.set_linewidth(.8); s.set_color("#444444")

    cb = fig.colorbar(im, ax=ax, orientation="vertical", pad=.012, fraction=.0265)
    cb.set_label("storms per 1$^{\\circ}$ cell, 1884 to 2023", fontsize=8.2)
    cb.ax.tick_params(labelsize=7.5, length=2.5, width=.6)
    cb.outline.set_linewidth(.7); cb.outline.set_edgecolor("#444444")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{dst}/Fig1_domain.{ext}", dpi=400, bbox_inches="tight")
    print("written:", dst)


if __name__ == "__main__":
    main()
