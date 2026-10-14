# TCCounts

Tropical cyclone count trends for the western North Pacific basin and the
Philippine Area of Responsibility (PAR), recomputed from IBTrACS under every
defensible combination of analytic choices, with the checks that reproduce every
reported value.

Four choices enter any such trend, and no convention settles them:

* **grading requirement**: which operational agencies must have graded a storm for
  it to count (JMA, JTWC, CMA, JMA or JTWC, any of the three, all three)
* **intensity threshold**: all systems, tropical depression or stronger, or
  tropical storm or stronger, each agency judged by its own categories
* **counting rule**: which track entries constitute a storm (4 options)
* **start year**: the first year of the fitted window (50 options, 1951 to 2000;
  every window ends in 2023)

6 x 3 x 4 x 50 = 3,600 specifications per domain, 7,200 for the basin and the PAR
together. HKO, any of the four agencies and the raw total are fitted as three
further series, 2,800 specifications, for 10,000 in all. Each slope is tested
against zero with a standard error corrected for serial dependence; the size of
that test is measured by simulation and calibrated against fitted noise models,
and the specification curve is also tested as a whole.

At tropical-storm strength the agencies largely agree, and the count falls in
both domains. Most of the decline in the total count lies in systems that no
agency rated as a tropical storm or that no agency graded, and those entries
changed in dated steps rather than gradually. PAGASA's own record of storms
entering the PAR is compared with the archive by `src/pagasa.py`.

## Reproducing it

Neither the IBTrACS archive nor PAGASA's records are redistributed. Get IBTrACS
v04r01 for the western North Pacific from NOAA NCEI and run everything from the
repository root.

### Version 3: the specification curve and its test

```
python src/clip_par.py   ibtracs.WP.list.v04r01.csv par_clipped.csv
python src/multiverse.py ibtracs.WP.list.v04r01.csv par_clipped.csv data/multiverse.csv
python src/coverage.py   ibtracs.WP.list.v04r01.csv par_clipped.csv data/coverage.csv
python src/map_domain.py ibtracs.WP.list.v04r01.csv figures/
python src/analysis.py   data/multiverse.csv par_clipped.csv ibtracs.WP.list.v04r01.csv
python src/validate.py
```

Steps 2 and 3 reproduce the committed CSVs. Across numpy and scipy versions the
no-trend bands in `multiverse.csv` can differ in the last digits (by about
1e-15); the slopes and verdicts do not change. Timings from a clean tree, single
core:

| script | what it does | time |
|---|---|---|
| `clip_par.py` | clips IBTrACS to the PAR hexagon | 9 s |
| `multiverse.py` | fits and tests the 2,400 specifications of version 3 | 31 s |
| `coverage.py` | annual counts and per-agency graded counts | 9 s |
| `map_domain.py` | Fig 1, the study area | 3 s |
| `analysis.py` | recomputes and checks every value of version 3 | 34 s |
| `validate.py` | size and power of the significance test, by simulation | 100 s |

`analysis.py` prints OK or MISMATCH beside each value and exits non-zero if any
check fails. A clean run is **199 checks, 0 failed**.

### Version 4: intensity threshold, storm classes, dating and calibrated tests

Of the scripts below, only `threshold.py` reads the archive; the others read the
per-storm flags it writes.

```
python src/threshold.py   ibtracs.WP.list.v04r01.csv par_clipped.csv data/
python src/decompose.py   data/storm_flags.csv data/
python src/dating.py      data/storm_flags.csv data/
python src/calibrate.py   data/storm_flags.csv data/
python src/curve_test.py  data/storm_flags.csv data/
python src/end_year.py    data/storm_flags.csv data/
python src/check_v4.py    data/
python src/fig_curves.py  data/grid_threshold.csv figures/
python src/fig_classes.py data/class_counts.csv figures/
python src/fig_design.py  figures/
```

With the library versions used for version 4, every committed CSV is reproduced
byte for byte; with other numpy versions the simulated outputs can differ, which
the checks allow for. Timings from a clean tree, single core, on the machine used
for version 4:

| script | what it does | time |
|---|---|---|
| `threshold.py` | per-storm flags; 10,000 specifications over requirement, threshold, counting rule and start year | 27 s |
| `decompose.py` | splits each trend by storm class; tropical-storm coverage by era; CMA wind sensitivity | 2 s |
| `dating.py` | step and smooth-trend models for coverage shares and class counts | 22 s |
| `calibrate.py` | test size under AR(1) and fitted AR(p) and ARMA(1,1) noise; power by window length; calibrated p-values | 25 s |
| `curve_test.py` | 2,000 block-bootstrap archives with no trend, run through every specification | 6 min |
| `end_year.py` | verdicts for end years 2005 to 2023 | 2 s |
| `check_v4.py` | checks every reported version 4 value | 4 s |
| `fig_curves.py` | Fig 2, fitted slope against start year | 5 s |
| `fig_classes.py` | Fig 3, annual counts by storm class with the dated steps | 13 s |
| `fig_design.py` | Online Resource 1, Fig S1, the enumeration | 2 s |

The three figure scripts were timed on a slower machine, on which `curve_test.py`
took 12.5 min.

`check_v4.py` works like `analysis.py`. A clean run is **376 checks, 0 failed**.
Deterministic values are checked to the rounding they were reported at. Simulated
values (test size, power, calibrated p-values, the curve-level test) are checked
with Monte Carlo tolerances, so the suite also passes where numpy's random streams
differ from the version used here.

### PAGASA's record

PAGASA's records are preliminary operational records, available from PAGASA on
request. They are not redistributed here, and neither is anything derived from
them. With PAGASA's list of tropical cyclones that crossed the PAR and its
six-hourly track file (described in the header of `src/pagasa.py`), write the
outputs to a directory outside the repository:

```
python src/pagasa.py       LIST.xlsx TRACKS.csv ibtracs.WP.list.v04r01.csv data/storm_flags.csv OUTDIR
python src/check_pagasa.py LIST.xlsx OUTDIR data/
python src/fig_classes.py  data/class_counts.csv OUTDIR OUTDIR/pagasa_annual.csv
```

`pagasa.py` runs in about 70 s and `check_pagasa.py` in about 27 s on the slower
machine. A clean run of `check_pagasa.py` is **101 checks, 0 failed**. Given
PAGASA's annual series, `fig_classes.py` adds it to the PAR panels, as in the
manuscript; the committed `figures/Fig3_classes.pdf` is drawn from IBTrACS alone.

## Requirements

`numpy`, `pandas`, `scipy` and `matplotlib`; `statsmodels` for `calibrate.py`
and for the two PAGASA scripts, which use it; `openpyxl` for the two PAGASA
scripts, which read a spreadsheet. No network access at run time: the shoreline the map
needs is committed in `data/`.

## Layout

```
src/     version 3: clip_par, multiverse, coverage, map_domain, analysis, validate
         version 4: classes, threshold, decompose, dating, calibrate, curve_test,
                    end_year, check_v4, pagasa, check_pagasa,
                    fig_curves, fig_classes, fig_design
data/    multiverse.csv, coverage.csv, coastline_wnp.npz, README.md
         version 4 outputs, listed in data/README.md
figures/ Fig1_domain.pdf, Fig2_curves.pdf, Fig3_classes.pdf, ESM_FigS1_design.pdf
```

`data/README.md` documents the columns, the PAR hexagon, and two column-name
traps in IBTrACS that will bite anyone who recomputes this from a shapefile
export. The figure scripts of version 3, `figures.py` and `schematic.py`, and
their figures are in the v3.0.0 release.

## License

MIT.
