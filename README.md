# BayesPrePost-FBST

R/C++ code for the tables and figures in `manuscripts/original/article.tex`,
including the appendix and Supporting Information.

## Run

Requires R 4.4 or later and a C++ compiler compatible with Rcpp
(Rtools on Windows). Install once:

```r
install.packages(c("Rcpp", "statmod"))
```

Run from the repository root:

```sh
Rscript main.R                          # all tables and figures
Rscript main.R tables                   # tables only
Rscript main.R figures                  # figures only
Rscript main.R vaping_data vaping_freq  # selected steps
```

In R, use `source("main.R")`. To run an individual function, first load
`source("R/load.R")`, then call, for example, `compute_power()`.

Tables are saved in `results/tables/` as unrounded CSVs and plain TeX
tabulars, named after the article labels: `tab:design` becomes `design.csv`
and `design.tex`. Figures are saved in `figures/` as PDF, SVG and PNG.
A full run can take many hours; `results/cache/` reuses completed calculations
and saves enumeration progress. `results/sessionInfo.txt` records the R environment.

## Structure

- `main.R`: entry point and selection of steps.
- `R/Data/`: aggregate vaping counts and [data provenance](R/Data/README.md).
- `R/Data_preparation/`: count loading and admissible reconstructions.
- `R/Method/`: priors, posterior, e-values, calibration and sampling distributions.
- `R/Tables/`: one script per table or related group.
- `R/Figures/`: figure scripts.
- `R/Utils/`: caching and export; `R/config.R`: input/output paths.

Study scenarios and numerical settings are in `R/Method/study_settings.R`.
Generated files, manuscripts and local working material are excluded from Git.

`main.R` maps each step to its calculation function. Scripts are named after
their calculation, such as `power.R`, `typeI_error.R` and `prior_sensitivity.R`;
output filenames follow the article's LaTeX labels.
Table 1 defines notation and needs no calculation.
`reconstruction` also writes `results/reconstruction_solutions.csv`.

| Step | Script in `R/Figures/` | Original figures |
|---|---|---|
| `diagrams` | `model_diagrams.R` | 1–2 |
| `intervals` | `posterior_intervals.R` | 3 |
| `errors` | `error_curves.R` | 4 |
| `densities` | `prior_posterior.R` | 5 |

Every step runs independently and reuses cached calculations. Figures 1–2
are redrawn schematics; the restriction diagram uses the informative prior.
