"""
Step 1. Clip the IBTrACS western North Pacific archive to the PAGASA
Philippine Area of Responsibility (PAR) hexagon.

Input   ibtracs.WP.list.v04r01.csv  (or the equivalent WP-basin CSV export)
        Obtain from the NOAA NCEI IBTrACS archive; see data/README.md.
Output  par_clipped.csv

The PAR is a six-vertex polygon, not a bounding box. The diagonal edge from
(15 N, 115 E) to (21 N, 120 E) matters: a bounding-box clip admits storms the
PAR excludes.

Usage   python src/clip_par.py ibtracs.WP.list.v04r01.csv par_clipped.csv
"""
import sys
import numpy as np
import pandas as pd
from matplotlib.path import Path

PAR_VERTICES = [(115, 5), (115, 15), (120, 21), (120, 25), (135, 25), (135, 5)]  # (lon, lat)


def clip(df, lon_col="LON", lat_col="LAT"):
    """Return the rows whose fix falls inside or on the PAR boundary."""
    poly = Path(PAR_VERTICES)
    lon = pd.to_numeric(df[lon_col], errors="coerce").values
    lat = pd.to_numeric(df[lat_col], errors="coerce").values
    ok = np.isfinite(lon) & np.isfinite(lat)
    inside = np.zeros(len(df), dtype=bool)
    # radius=-1e-9 dilates the polygon by a hair so fixes reported exactly on a
    # boundary meridian or parallel are kept. IBTrACS reports to 0.1 degrees, so
    # a strict test silently drops well over a thousand legitimate fixes.
    inside[ok] = poly.contains_points(np.c_[lon[ok], lat[ok]], radius=-1e-9)
    return df[inside]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "ibtracs.WP.list.v04r01.csv"
    dst = sys.argv[2] if len(sys.argv) > 2 else "par_clipped.csv"
    df = pd.read_csv(src, low_memory=False, skiprows=[1])  # row 1 is the units row
    out = clip(df)
    out.to_csv(dst, index=False)
    print(f"{len(df):,} fixes in -> {len(out):,} fixes inside the PAR")
    print(f"{out.SID.nunique():,} unique storms, seasons {int(out.SEASON.min())}-{int(out.SEASON.max())}")
    print(f"written: {dst}")


if __name__ == "__main__":
    main()
