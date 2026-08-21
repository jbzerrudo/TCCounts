# Data

## What is here

| file | what it is |
|---|---|
| `multiverse.csv` | one row per fitted specification, 2,400 rows |
| `coverage.csv` | annual main-track storm count and per-agency graded count, both domains, 1884 to 2023 |
| `coastline_wnp.npz` | GSHHS shoreline for the map frame, so the figure script needs no network |

### `multiverse.csv`

| column | meaning |
|---|---|
| `domain` | `WNP basin` or `PAR` |
| `rule` | counting rule: `all SIDs`, `main track`, `>=2 fixes`, `tropical nature` |
| `agency` | grading requirement: `raw total`, `JTWC`, `JMA`, `JMA+JTWC`, `JMA+JTWC+CMA`, `all four` |
| `start` | first year of the fitted window; every window ends in 2023 |
| `slope` | OLS slope, storms per year |
| `lo`, `hi` | the two-sided 95% no-trend band from the serial-dependence corrected standard error |
| `verdict` | `-` significant decline, `0` no detectable trend, `+` significant increase |

6 requirements x 4 rules x 50 start years x 2 domains = 2,400 rows exactly.

### `coverage.csv`

`SEASON`, `domain`, `total`, then `JMA`, `JTWC`, `CMA`, `HKO`. `total` is every
main-track storm in that season and domain; each agency column is how many of
those storms that agency assigned an intensity to. An agency column can exceed
none of them and never exceeds `total`. The ratio is the coverage share whose
growth the paper is about.

### `coastline_wnp.npz`

Shoreline polylines for 99 to 181 E and -1 to 56 N, held as a flat `pts` array
of float32 lon/lat pairs plus a `lens` index giving each segment's length.
Extracted once from the GSHHS shoreline distributed with `basemap-data`, so
`src/map_domain.py` runs offline with only matplotlib, numpy and pandas.

## What is not here, and why

The clipped best-track archive is roughly 35 MB and is not redistributed. It is
derived from IBTrACS, which NOAA maintains and which should be obtained from the
source so you get the current revision.

1. Download the **western North Pacific (WP) basin, version 4 revision 1** CSV from
   the NOAA NCEI IBTrACS archive. Cite Knapp et al. (2010),
   <https://doi.org/10.1175/2009BAMS2755.1>. Confirm the current filename at the
   archive; the naming convention is `ibtracs.WP.list.v04r01.csv`.
2. `python src/clip_par.py ibtracs.WP.list.v04r01.csv par_clipped.csv`
3. `python src/multiverse.py ibtracs.WP.list.v04r01.csv par_clipped.csv data/multiverse.csv`
4. `python src/coverage.py ibtracs.WP.list.v04r01.csv par_clipped.csv data/coverage.csv`

Steps 3 and 4 should reproduce `multiverse.csv` and `coverage.csv` byte for byte.
If they do not, the IBTrACS revision has changed since this was run, which is
itself worth knowing.

## The PAR is a hexagon

Vertices, from PAGASA: (5 N, 115 E), (15 N, 115 E), (21 N, 120 E), (25 N, 120 E),
(25 N, 135 E), (5 N, 135 E). The diagonal from (15 N, 115 E) to (21 N, 120 E) is
load-bearing. A bounding-box clip at 5-25 N, 115-135 E admits 129 main-track
storms the PAR excludes, 4.7% of the 2,757 the polygon contains, and the two
clips have identical coordinate extents, so bounds alone will not tell you which
one you have. Test a corner: (25 N, 115 E) is inside the box and outside the PAR.

One further trap. IBTrACS reports positions to 0.1 degrees, so many fixes land
exactly on a boundary meridian or parallel. A strict point-in-polygon test drops
them. `clip_par.py` dilates the polygon by 1e-9 degrees to keep them.

## Two column-name traps

The NCEI CSV names the JMA grade column `TOKYO_GRADE`. Shapefile exports truncate
it to `TOK_GRADE`. Every script here accepts either.

IBTrACS `NATURE` uses `TS` for tropical and `DS` for disturbance. `DS` is not a
tropical stage, so the tropical-nature counting rule tests for `TS` alone.
