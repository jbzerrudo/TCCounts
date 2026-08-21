# TCCounts

Every defensible specification of a tropical cyclone count trend in the western
North Pacific basin and the Philippine Area of Responsibility, enumerated and
fitted, and the reason they disagree.

Three choices enter any such trend and no convention settles any of them:

* **grading requirement** which operational agencies must have graded a storm for
  it to count (6 options, including the raw total that asks for no grade at all)
* **counting rule** which track entries constitute a storm (4 options)
* **start year** the first year of the fitted window (50 options, 1951 to 2000;
  every window ends in 2023)

6 x 4 x 50 = 1,200 specifications per domain, fitted for the basin and again for
the PAR polygon, so 2,400 in total. Each slope is tested against zero using a
standard error corrected for serial dependence, whose size and power are measured
by simulation rather than assumed.

The archived count falls in both domains. The specifications that return no trend
are the ones restricted to an agency whose coverage of the archive was still
growing, because a count restricted to such an agency gains entries from
bookkeeping and rises against a falling total.

## Reproducing it

The clipped archive is not redistributed. Get IBTrACS v04r01 for the western
North Pacific from NOAA NCEI, then:

```
python src/clip_par.py   ibtracs.WP.list.v04r01.csv par_clipped.csv
python src/multiverse.py ibtracs.WP.list.v04r01.csv par_clipped.csv data/multiverse.csv
python src/coverage.py   ibtracs.WP.list.v04r01.csv par_clipped.csv data/coverage.csv
python src/map_domain.py ibtracs.WP.list.v04r01.csv figures/
python src/figures.py    par_clipped.csv ibtracs.WP.list.v04r01.csv data/multiverse.csv figures/
python src/analysis.py   data/multiverse.csv par_clipped.csv ibtracs.WP.list.v04r01.csv
python src/validate.py
```

Steps 2 and 3 reproduce the committed CSVs byte for byte. Timings from a clean
tree, single core:

| script | what it does | time |
|---|---|---|
| `clip_par.py` | clips IBTrACS to the PAR hexagon | 9 s |
| `multiverse.py` | fits and tests all 2,400 specifications | 31 s |
| `coverage.py` | annual counts and per-agency graded counts | 9 s |
| `map_domain.py` | Fig 1, the study area | 3 s |
| `figures.py` | Fig 2 counts and coverage, Fig 3 specification curve | 11 s |
| `analysis.py` | recomputes and checks every value in the manuscript | 34 s |
| `validate.py` | size and power of the significance test, by simulation | 100 s |

`analysis.py` prints OK or MISMATCH beside each value and exits non-zero if any
check fails. A clean run is **193 checks, 0 failed**.

## Requirements

`numpy`, `pandas`, `scipy`, `matplotlib`. Nothing else, and no network access at
run time: the shoreline the map needs is committed in `data/`.

## Layout

```
src/    clip_par, multiverse, coverage, map_domain, figures, analysis, validate
data/   multiverse.csv, coverage.csv, coastline_wnp.npz, README.md
figures/ Fig1_domain.pdf, Fig1.pdf, Fig2.pdf
```

`data/README.md` documents the columns, the PAR hexagon, and two column-name
traps in IBTrACS that will bite anyone who recomputes this from a shapefile
export.

## Licence

MIT.
