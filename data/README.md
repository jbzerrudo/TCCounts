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
none of them and never exceeds `total`. The ratio is that agency's coverage share,
whose changes `dating.py` dates.

### `coastline_wnp.npz`

Shoreline polylines for 99 to 181 E and -1 to 56 N, held as a flat `pts` array
of float32 lon/lat pairs plus a `lens` index giving each segment's length.
Extracted once from the GSHHS shoreline distributed with `basemap-data`, so
`src/map_domain.py` runs offline with only matplotlib, numpy and pandas.

## Version 4 outputs

Thresholds, each agency judged by its own categories at any fix inside the domain:
`all` is the definition in `multiverse.csv` (JMA grade 2-5 or 9, or any positive
wind from JTWC, CMA or HKO); `TD` is tropical depression or stronger (JMA grade
2-5 or 9, `USA_SSHS` >= -1, `CMA_CAT` 1-6, `HKO_CAT` TD or stronger); `TS` is
tropical storm or stronger (JMA grade 3-5 or 9, `USA_SSHS` >= 0, `CMA_CAT` 2-6,
`HKO_CAT` TS or stronger). Requirements: `JMA`, `JTWC`, `CMA`, `JMA|JTWC`,
`any of 3`, `all 3` (the three basin-wide agencies), with `HKO`, `any of 4` and
`raw total` kept as separate series.

| file | what it is |
|---|---|
| `storm_flags.csv` | one row per storm, domain and counting rule: season and a 0/1 flag for each agency and threshold, plus `CMA_TSwind` (CMA by wind, >= 34 kt) |
| `grid_threshold.csv` | one row per specification, 10,000 rows: `domain`, `rule`, `requirement`, `threshold`, `start`, `slope`, `lo`, `hi`, `verdict` as in `multiverse.csv` |
| `archive_meta.csv` | `TRACK_TYPE` row counts for 1884-2023 and after 2023 |
| `cma_consistency.csv` | CMA fixes whose category disagrees with their wind |
| `class_counts.csv` | annual main-track counts per storm class, both domains, 1884-2023 |
| `decomposition.csv` | slope, verdict and share of the total slope per class, windows from 1951, 1961, 1977 |
| `era_means.csv` | storms per year per class in 1951-1976, 1977-2000, 2001-2023 |
| `ts_coverage.csv` | share of tropical-storm-strength storms each agency also rates tropical storm or stronger, by era |
| `cma_wind_sensitivity.csv` | CMA classified by wind instead of category |
| `dating_shares.csv`, `dating_counts.csv` | selected model, break years, support intervals and segment means for coverage shares and class counts |
| `size_ar1.csv`, `size_fitted.csv` | rejection rate with no trend under AR(1), and under the AR(p) chosen by AIC and an ARMA(1,1) |
| `power.csv` | trend reaching 80% power, and the change it implies as a share of the window mean, windows of 24-73 years |
| `calibrated_p.csv` | two-sided p-values of observed slopes under the fitted noise models |
| `curve_test.csv` | observed significant declines per family against 2,000 no-trend archives, block lengths 5, 10, 15 |
| `end_year.csv` | share of windows with a significant decline, and the slope from 1951, for end years 2005-2023 |

Storm classes in `class_counts.csv`: `ts_or_stronger`, `depression_only`,
`below_depression` (a positive wind but no depression category from any agency)
and `no_intensity` partition `total`; `jma_graded`, `jtwc_not_jma`, `cma_hko_only`
and `hko_only` split it by which agencies graded the storm; `ts_all3`, `ts_some3`
and `ts_hko_only` split `ts_or_stronger` by agreement among JMA, JTWC and CMA.

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

Steps 3 and 4 should reproduce `multiverse.csv` and `coverage.csv`; across numpy
and scipy versions the `lo` and `hi` bands can differ in the last digits. If the
slopes or verdicts differ, the IBTrACS revision has changed since this was run,
which is itself worth knowing.

PAGASA's tropical cyclone records are not here either, and neither is anything
`src/pagasa.py` derives from them. They are preliminary operational records,
available from PAGASA on request; the main README shows how to run the comparison
with the outputs kept outside the repository.

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
