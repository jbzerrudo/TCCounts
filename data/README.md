# Data

## What is here

`multiverse.csv` — every fitted specification, 2,400 rows.

| column | meaning |
|---|---|
| `domain` | `WNP basin` or `PAR` |
| `rule` | counting rule: `all SIDs`, `main track`, `>=2 fixes`, `tropical nature` |
| `agency` | grading requirement: `raw total`, `JTWC`, `JMA`, `JMA+JTWC`, `JMA+JTWC+CMA`, `all four` |
| `start` | first year of the fitted window; every window ends in 2023 |
| `slope` | ordinary least squares slope, storms per year |
| `lo`, `hi` | 95% no-trend band from the effective-sample-size correction |
| `verdict` | `-` significant decline, `+` significant increase, `0` no detectable trend |

`coverage.csv` — annual counts and per-agency coverage, 1884 to 2023, both domains,
main-track storms only. This is what Table 1 and Figure 1 are built from.

| column | meaning |
|---|---|
| `SEASON` | year |
| `domain` | `WNP basin` or `PAR` |
| `total` | unique main-track storms with at least one fix inside the domain |
| `JMA` | of those, the number JMA graded (TOKYO_GRADE in 2, 3, 4, 5, 9) |
| `JTWC` | of those, the number with a positive USA_WIND |
| `CMA` | of those, the number with a positive CMA_WIND |
| `HKO` | of those, the number with a positive HKO_WIND |

Together these two files reproduce both tables, both figures, and every percentage
quoted in the manuscript, without redownloading the archive.

## What is not here, and why

The clipped best-track archive, `par_clipped.csv`, is about 34 MB and is not
redistributed. It is derived from IBTrACS, which NOAA maintains and which should be
obtained from the source so you get the current revision.

1. Download the **western North Pacific (WP) basin, version 4 revision 1** CSV from
   <https://www.ncei.noaa.gov/products/international-best-track-archive>
2. Set `WNP` to its path, then from the repository root:

```
python src/clip_par.py   "$WNP" par_clipped.csv
python src/multiverse.py "$WNP" par_clipped.csv data/multiverse.csv
python src/figures.py    par_clipped.csv "$WNP" data/multiverse.csv figures/
python src/analysis.py   data/multiverse.csv par_clipped.csv "$WNP"
```

The last line prints `137 checks run, 0 failed`. Windows `cmd` equivalents are in the
repository README.

`src/validate.py` establishes the size and power of the significance test by
simulation. It takes a few minutes and is not part of the 137 checks.

## Superseded

Earlier versions of this repository carried `par_annual_series.csv`,
`start_year_sweep.csv`, `backtest.json` and a classified-versus-unclassified
decomposition. That analysis was replaced at v2.0.0. The files remain available in the
v1.4.0 release and its Zenodo record.
