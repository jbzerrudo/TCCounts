# TCCounts

Annual counts of tropical cyclones entering the Philippine Area of Responsibility
(PAR) are reported to be declining. This repository shows that the direction and
significance of that trend are set by the analyst's choice of start year acting on
an archive of changing completeness, and gives you everything needed to check the
claim from the public IBTrACS record.

## The result in one table

Split each year's PAR storms into two classes: those an operational agency ever
classified for intensity, and those carried as a position on a map and nothing else.
The two classes partition the total, so their trends sum to the total's trend exactly.

Fit the raw count over windows that all end in 2023 but begin anywhere from 1884 to 2000:

| | range of OLS slope | significant increase | significant decrease |
|---|---|---|---|
| **all track entries** | +0.058 to −0.146 storms yr⁻¹ | 4 of 59 start years | 16 of 59 start years |
| **classified only** | −0.056 to +0.101 storms yr⁻¹ | 0 of 25 | 0 of 25 |

Same archive, same domain, same estimator, opposite conclusions.

The mechanism is compositional:

| period | all entries | classified | unclassified | unclassified share |
|---|---|---|---|---|
| 1884–1922 | 14.9 | 0.0 | 14.9 | 100% |
| 1923–1944 | 19.3 | 0.0 | 19.3 | 100% |
| 1945–1950 | 21.0 | 14.7 | 6.3 | 30% |
| 1951–1976 | 25.1 | 18.6 | 6.5 | 26% |
| 1977–2000 | 22.1 | 20.0 | 2.1 | 10% |
| 2001–2023 | 19.3 | 18.9 | 0.4 | 2% |

No PAR storm before 1945 carries intensity information from any agency: zero of 1,004.

Spur tracks are excluded. IBTrACS gives a diverging agency track its own storm identifier, so
counting every SID double-counts those storms; 43 PAR identifiers are spur-only.
JTWC reporting begins in 1945, JMA grading in 1951, WMO winds in 1977.

## Reproducing it

```bash
pip install -r requirements.txt

# 1. get IBTrACS WP v04r01 from NOAA NCEI (see data/README.md), then:
python src/clip_par.py ibtracs.WP.list.v04r01.csv par_clipped.csv
python src/build_series.py par_clipped.csv data/par_annual_series.csv

# 2. recompute and check every number quoted in the manuscript
python src/analysis.py data/par_annual_series.csv data/start_year_sweep.csv

# 3. redraw the figure
python src/figure.py data/par_annual_series.csv figures/
```

`analysis.py` prints `OK` or `MISMATCH` beside each value and runs 58 checks. A clean
run prints no `MISMATCH`. If you get one, the IBTrACS revision has moved and the
manuscript numbers need restating; shipping the check rather than the result is the point.

Step 1 needs the raw archive. Steps 2 and 3 run from `data/` alone, which is in this
repository, so the analysis is verifiable without downloading anything.

## Layout

```
src/clip_par.py        IBTrACS WP  ->  PAR-clipped fixes (hexagon, not bounding box)
src/build_series.py    PAR-clipped ->  annual series, 1884-2023
src/analysis.py        every manuscript number, each checked
src/figure.py          Figure 1
data/                  derived series, start-year sweep, backtest, provenance notes
figures/               Figure 1 as PDF
```

## Citing

See `CITATION.cff`. The archived release carries a DOI via Zenodo.

## License

MIT for the code. The derived series in `data/` are a transformation of IBTrACS,
which is a US Government work in the public domain; cite Knapp et al. (2010).
