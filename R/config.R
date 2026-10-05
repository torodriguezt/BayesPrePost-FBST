# Input and output paths. Run every script from the repository root.
# Priors, scenarios and numerical controls are in Method/priors.R.

FBST_PATHS <- list(
  vaping_data = "R/Data/vaping_mccauley2023.csv",
  results = "results",
  figures = "figures"
)

fbst_output_dir <- function() {
  dir.create(FBST_PATHS$results, recursive = TRUE, showWarnings = FALSE)
  FBST_PATHS$results
}

fbst_figure_dir <- function() {
  dir.create(FBST_PATHS$figures, recursive = TRUE, showWarnings = FALSE)
  FBST_PATHS$figures
}

fbst_vaping_path <- function() FBST_PATHS$vaping_data
