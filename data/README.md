# Data

## What is here

`par_annual_series.csv` — the derived annual series, 1884 to 2023, one row per season.
`start_year_sweep.csv` — OLS slope and block-bootstrap no-trend band for every start year.
`backtest.json` — rolling-origin forecast comparison.

| column | meaning |
|---|---|
| `SEASON` | year |
| `total` | every unique storm with at least one fix inside the PAR |
| `any_int` | **classified**: a JMA grade **or** a JTWC wind at any PAR fix |
| `pos_only` | **unclassified**: `total - any_int`; a position on a map and nothing else |
| `jma_class` | the stricter JMA-graded-only subset, for sensitivity |

`total = any_int + pos_only` exactly, every year. That identity is the whole argument.

## What is not here, and why

The clipped best-track archive is roughly 43 MB and is not redistributed. It is
derived from IBTrACS, which NOAA maintains and which should be obtained from the
source so you get the current revision.

1. Download the **western North Pacific (WP) basin, version 4 revision 1** CSV from
   the NOAA NCEI IBTrACS archive. Cite Knapp et al. (2010),
   <https://doi.org/10.1175/2009BAMS2755.1>. Confirm the current filename at the
   archive; the naming convention is `ibtracs.WP.list.v04r01.csv`.
2. `python src/clip_par.py ibtracs.WP.list.v04r01.csv par_clipped.csv`
3. `python src/build_series.py par_clipped.csv data/par_annual_series.csv`

Step 2 should reproduce `par_annual_series.csv` byte for byte. If it does not, the
IBTrACS revision has changed since this was run, which is itself worth knowing.

## The PAR is a hexagon

Vertices, from PAGASA: (5 N, 115 E), (15 N, 115 E), (21 N, 120 E), (25 N, 120 E),
(25 N, 135 E), (5 N, 135 E). The diagonal from (15 N, 115 E) to (21 N, 120 E) is
load-bearing. A bounding-box clip at 5-25 N, 115-135 E admits storms the PAR
excludes, and the two clips have identical coordinate extents, so bounds alone will
not tell you which one you have. Test a corner: (25 N, 115 E) is inside the box and
outside the PAR.

One further trap. IBTrACS reports positions to 0.1 degrees, so many fixes land
exactly on a boundary meridian or parallel. A strict point-in-polygon test drops
them. `clip_par.py` dilates the polygon by 1e-9 degrees to keep them.
