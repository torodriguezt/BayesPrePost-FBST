# Bayesian pretest–posttest analysis for binary outcomes

Code for the paper *"Bayesian analysis for pretest-posttest binary outcomes
with adaptive significance levels"*. The method compares two dependent
proportions (pre- vs. post-treatment success probabilities) with an
Olkin–Liu bivariate beta prior, which yields a closed-form bivariate
beta-binomial posterior. Precise hypotheses are tested with the Full Bayesian
Significance Test (FBST); the significance threshold `k*` is adaptive
(a function of the sample size and the data) and the loss weight `a*` is
calibrated prior-based so that the Type-I Bayes risk is approximately 0.05.

## Requirements

R with packages: `Rcpp`, `ALA` (TVSFP data), `dplyr`, `tidyr`, `ggplot2`,
`lattice`, plus a C++ compiler for `Rcpp::sourceCpp`.

## Reproducing the paper

From the repository root:

```r
Rscript run_all.R
```

This regenerates every table and figure (about 2 hours). Outputs are written
to `Figures/` (PNG) and `output/` (CSV / RDS / TeX fragments); both are
git-ignored. Each step can also be run on its own.

## Structure

| Path | Purpose |
|------|---------|
| `run_all.R` | Master script: runs the full pipeline in order |
| `src/BivBetaBinom.cpp` | Computational core: closed-form posterior, FBST e-values by 2D quadrature, adaptive cutoff `k*`, prior- and posterior-based samplers |
| `R/01_priors_config.R` | Hyperparameters of the three priors (single source of truth) |
| `R/02_fit_priors_kl.R` | Derives the priors by Kullback–Leibler minimisation |
| `R/03_calibrate_weights.R` | Prior-based calibration of the loss weight `a*` as a function of `n` and the prior (reference lookup table) |
| `R/04_simulation_estimation_figures.R` | Estimation study: prior/posterior surfaces, posterior mean and mode vs. `n` |
| `R/05_simulation_prior_sensitivity.R` | Prior-sensitivity of the FBST decision on simulated data (heatmaps by prior category) |
| `R/06_application_tvsfp.R` | TVSFP application: FBST tables under the three priors (calibrated `a*`) and McNemar benchmark |
| `R/07_application_figures.R` | TVSFP posterior surfaces and error curves |

## Data

The application uses the Television School and Family Smoking Prevention and
Cessation Project (TVSFP) subset shipped with the R package `ALA`
(4 Los Angeles schools; binary outcome: THKS score >= 3). No data files are
stored in this repository.
