# TCCounts

Annual counts of tropical cyclones entering the Philippine Area of Responsibility
(PAR) show a decline since 1951. This repository shows that the decline is a
property of the archive rather than of the climate, and gives you everything
needed to check that claim from the public IBTrACS record.

## The result in one table

Split each year's PAR storms into two piles: those an operational agency graded for
intensity, and those carried as a position on a map and nothing else.

| period | all entries | intensity-classified | unclassified | unclassified share |
|---|---|---|---|---|
| 1951-1976 | 25.8 | 12.4 | 7.0 | 27% |
| 1977-2000 | 22.8 | 17.2 | 2.4 | 10% |
| 2001-2023 | 19.3 | 16.2 | 0.4 | 2% |

Over 1951 to 2023 the raw count falls at 0.122 storms per year, which a moving-block
bootstrap places outside the no-trend null. The unclassified component alone falls at
0.132 per year, more than accounting for the whole decline. The classified count does
not fall in any window tested.

Over the full 1923 to 2023 record the raw count looks flat, at -0.0070 storms per
year. That flatness is a cancellation, not a null: the classified component rises at
+0.2082 and the unclassified falls at -0.2152, and the two sum to the observed value
by construction. Both components are among the most significant trends in the data
and neither is physical.

No PAR fix before 1945 carries intensity information from any agency. JTWC reporting
begins in 1945, JMA grading in 1951, WMO winds in 1977.

## Reproducing it

```bash
pip install -r requirements.txt

# 1. get IBTrACS WP v04r01 from NOAA NCEI (see data/README.md), then:
python src/clip_par.py ibtracs.WP.list.v04r01.csv par_clipped.csv
python src/build_series.py par_clipped.csv data/par_annual_series.csv

# 2. recompute and check every number quoted in the manuscript
python src/analysis.py data/par_annual_series.csv par_clipped.csv

# 3. redraw the figure
python src/figure.py data/par_annual_series.csv figures/
```

`analysis.py` prints `OK` or `MISMATCH` beside each value. A clean run prints no
`MISMATCH`. If you get one, the IBTrACS revision has moved and the manuscript
numbers need restating; that is the point of shipping the check rather than the
result.

Step 1 needs the raw archive. Steps 2 and 3 run from `data/par_annual_series.csv`
alone, which is in this repository, so you can verify the analysis without
downloading anything.

## Layout

```
src/clip_par.py       IBTrACS WP  ->  PAR-clipped fixes
src/build_series.py   PAR-clipped ->  three annual series
src/analysis.py       every manuscript number, each checked
src/figure.py         Figure 1
data/                 the derived series, plus provenance notes
figures/              Figure 1 as PDF
```

## Citing

See `CITATION.cff`. The archived release carries a DOI via Zenodo.

## License

MIT for the code. The derived series in `data/` are a transformation of IBTrACS,
which is a US Government work in the public domain; cite Knapp et al. (2010).
