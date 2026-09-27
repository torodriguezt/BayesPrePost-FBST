# Adaptive FBST for pretest–posttest proportions

Code to reproduce the tables and figures of the article. It tests H: θ₁ = θ₂
with the FBST, a bivariate beta prior (Olkin–Liu) and an adaptive cutoff k*
that minimises α + β, computed by exact enumeration of the sample space.

## Installation

```r
install.packages(c("Rcpp", "statmod", "ALA"))
```

A C++ compiler for Rcpp is required.

## Usage

From the repository root:

```sh
Rscript run_all.R --stage=all --profile=full --keep-going
```

- Stages: `validate`, `pilot`, `main`, `boundaries`, `applications`, `vaping`, `appendix`, `export`.
- A single task: `--stage=main --task=power`.
- `--profile=pilot` runs a reduced test version.
- If interrupted, rerun the same command: cached results are reused.

The `vaping` stage analyses `vapeo_mccauley2023.csv` (McCauley et al., 2023):

```sh
Rscript run_all.R --stage=vaping --profile=full
```

## Outputs

- Tables: `output/reproducible_r/2.0.0/full/`
- Figures: `Figures/reproducible_r/2.0.0/full/`

## Code

| File | Contents |
|---|---|
| `R/exact_fbst.R`, `R/exact_fbst_core.cpp` | Posterior, e-values, predictive distributions and cutoffs |
| `R/01_priors_config.R`, `R/02_fit_priors_kl.R` | Priors and the KL-optimal prior |
| `R/pipeline_helpers.R` | Caching, calibration and export |
| `R/validate_numerics.R` | Numerical validation |
| `R/manuscript_experiments.R` | Numerical evaluation and appendix |
| `R/manuscript_applications.R` | TVSFP and onychomycosis applications |
| `R/application_vaping.R` | Vaping application |
| `R/04_simulation_estimation_figures.R` | Estimation appendix |
| `R/reproduction_report.R` | Run report |
