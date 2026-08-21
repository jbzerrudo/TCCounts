# TCCounts

Specification uncertainty in western North Pacific and Philippine tropical cyclone
count trends.

Fitting a linear trend to an annual tropical cyclone count requires three choices
that no convention settles:

| choice | options |
|---|---|
| **agency set** — which agencies must have graded a storm for it to count | JTWC; JMA; JMA+JTWC; +CMA; all four; raw total |
| **counting rule** — which track entries constitute a storm | all SIDs; main track; main track with 2 or more fixes; main track, tropical nature |
| **start year** — the first year of the fitted window | 1951 to 2000, all windows ending 2023 |

6 x 4 x 50 = 1,200 specifications per domain, fitted for the western North Pacific
basin and again for the PAGASA Philippine Area of Responsibility (PAR), so 2,400
in total. Each is tested against its own circular block bootstrap null.

**Result.** 42% of basin specifications and 30% of PAR specifications return a
significant decline; the rest return no detectable trend. Not one of the 2,400
returns a significant increase. Agency choice dominates the outcome; the counting
rule barely matters.

## Pipeline

Set `WNP` to wherever you saved the IBTrACS archive; it is not in this repository
and does not have to sit in the repository root.

Linux or macOS:

```
WNP=/path/to/ibtracs.WP.list.v04r01.csv

python src/clip_par.py   "$WNP" par_clipped.csv
python src/multiverse.py "$WNP" par_clipped.csv data/multiverse.csv
python src/figures.py    par_clipped.csv "$WNP" data/multiverse.csv figures/
python src/analysis.py   data/multiverse.csv par_clipped.csv "$WNP"
```

Windows `cmd`:

```
set WNP=C:\path\to\ibtracs.WP.list.v04r01.csv

python src\clip_par.py   "%WNP%" par_clipped.csv
python src\multiverse.py "%WNP%" par_clipped.csv data\multiverse.csv
python src\figures.py    par_clipped.csv "%WNP%" data\multiverse.csv figures/
python src\analysis.py   data\multiverse.csv par_clipped.csv "%WNP%"
```

The output directory argument to `figures.py` must end in a separator, and forward
slashes are safest on Windows, since a trailing backslash inside quotes escapes the
quote in `cmd`.

The last line must print `96 checks run, 0 failed`. Every number quoted in the
manuscript is recomputed there and compared against the value in the text.

`multiverse.py` takes about 40 seconds. The whole pipeline runs in under two minutes.

## Input

`ibtracs.WP.list.v04r01.csv`, the western North Pacific CSV export of IBTrACS
version 4 revision 1, from NOAA NCEI:
<https://www.ncei.noaa.gov/products/international-best-track-archive>

It is not redistributed here, and `par_clipped.csv` is a derived intermediate that
`.gitignore` excludes. Both regenerate from the command above.

## Files

```
src/clip_par.py    clips IBTrACS to the PAR hexagon (a six-vertex polygon,
                   not a bounding box; the boundary is treated as inside)
src/multiverse.py  enumerates and fits all 2,400 specifications
src/figures.py     regenerates Fig 1 and Fig 2
src/analysis.py    recomputes every manuscript value and checks it

data/multiverse.csv  one row per specification: domain, rule, agency, start
                     year, slope, bootstrap band, verdict
figures/Fig1.pdf     annual counts and agency coverage, both domains
figures/Fig2.pdf     slope against start year, one panel per agency set
```

## Notes on counting

IBTrACS assigns a separate identifier to an agency track that diverges from the
primary one, and its documentation lists such spurs among the reasons to exercise
care when counting. That is why the counting rule is enumerated rather than fixed.

JMA grades 2, 3, 4, 5 and 9 denote Tropical Depression, Tropical Storm, Severe
Tropical Storm, Typhoon and a tropical cyclone of Tropical Storm intensity or
higher. Grade 6 is extratropical and grade 7 marks a system just entering the JMA
area of responsibility, so neither counts as a grade here.

The NCEI CSV names the column `TOKYO_GRADE`; some shapefile exports shorten it to
`TOK_GRADE`. Both are accepted.

## Version history

- **v2.0.0** — specification-curve analysis. Supersedes v1.
- **v1.4.0** — classified-versus-unclassified decomposition of the PAR count.
  Preserved on Zenodo; its scripts were removed at v2.0.0.

## Citation

Archived at <https://doi.org/10.5281/zenodo.22026010> (concept DOI, always resolves
to the latest version). See `CITATION.cff`.

## License

MIT.
